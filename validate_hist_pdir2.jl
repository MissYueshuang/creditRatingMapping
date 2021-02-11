## Developed: 20200506
## if input is yyyymm, then is Historical; if input is yyyymmdd, then is revision
using Distributed # addprocs(6)
using SharedArrays
using CSV

mpath = @__DIR__
cd(mpath)
k = findfirst("VT\\", mpath)
CommonToolPath = string(mpath[1:k[end]], "VT_common_tools\\Julia_common_tools\\")
include(string(CommonToolPath,"Include_All_Functions.jl"))
Include_All_Functions(CommonToolPath)
GC_VT = global_constants_VT()
include(string(GC_VT["vtEnv"],"CRI_VT_MonitoringPDiR\\PDiR2.0\\calculate_hist_pdir2.jl"))

function validate_hist_pdir2(yyyymm,Econ=[])
    if isempty(Econ)
        econList = getAllEcons()
    else
        econList = Econ
    end
#    result = Dict()

    n = length(econList)
    failedEcon = SharedArray{Float64, 1}(n)

    @sync @distributed for econ = econList
        iEcon = indexin(econ,econList)[1]
#        println("This is econ $iEcon")
        pmtpath = string(GC_VT["dataEnv"], "ProductionData\\Historical\\$yyyymm\\Daily\\Products\\P17_Pdir2.0\\resultRating_$econ.mat")
        pdpath = string(GC_VT["dataEnv"], "ProductionData\\Historical\\$yyyymm\\Daily\\Products\\P2_Pd\\resultPD_($econ)_12.mat")
        if !isfile(pmtpath)
            if !isfile(pdpath)
#                println("There is no pdir2.0 file for econ $econ")
                continue
            else
                println("pd and pdir2.0 does not match for econ $econ")
            end
        end
        try
            pmt_pdir = vt_read_mat(pmtpath)
            if isempty(pmt_pdir)
                println("there is no data in file econ $econ")
                continue
            end
            vt_pdir = calculate_hist_pdir2(econ,yyyymm)
            out = value_Comparison(pmt_pdir,vt_pdir)
            if isempty(out)
                println("validation passed for econ $econ")
            else
                println("validation failed for econ $econ !!!")
                failedEcon[iEcon] = econ
    #            result[econ] = out
            end
        catch e
            println("Wrong with econ $econ")
        end
    end

    return failedEcon
end
