## Developed: 20200506
using CSV

mpath = @__DIR__
cd(mpath)
k = findfirst("VT\\", mpath)
CommonToolPath = string(mpath[1:k[end]], "VT_common_tools\\Julia_common_tools\\")
include(string(CommonToolPath,"Include_All_Functions.jl"))
Include_All_Functions(CommonToolPath)
GC_VT = global_constants_VT()
include(string(GC_VT["vtEnv"],"CRI_VT_MonitoringPDiR\\PDiR2.0\\calculate_daily_pdir2.jl"))

## This function is to validate incremental pdir2.0 since 20200101 to till_date
function validate_incre_PDiR2(till_date, Econ=[])
    path = string(GC_VT["dataEnv"], "ProductionData\\Recent\\Daily\\")
    folder = readdir(path)
    if isempty(Econ)
        econList = getAllEcons()
    else
        econList = Econ
    end
    result = Dict()
    # n = length(econList)
    # folder_new = folder[(startswith.(folder,"202009")) .| (startswith.(folder,"202010"))]
    folder_new = folder[occursin.(string(Int(floor(firstDailyDate/10000))),folder)]
    sort!(folder_new)
    idx = findall(startswith.(folder_new,string(firstDailyDate)))[1]
    folder_new = folder_new[idx:end]
    # till_folder_name = folder[occursin.(string(till_date), folder)][1]
    # m = findall(folder_new.==till_folder_name)[1]
    # failedEcon = SharedArray{Float64, 2}((n,m))

    for econ = econList
        result[econ] = Dict()
        for date = folder_new
            yyyymmdd = parse(Int64,date[1:8])
#            println("This is $yyyymmdd")
            pmtpath = string(path,date,"\\Products\\P17_Pdir2.0\\pdirNum_$econ.mat")
            if !isfile(pmtpath)
#                println("There is no file in Folder for $econ on date $yyyymmdd")
                continue
            end
            pmtdata = vt_read_mat(pmtpath)
            pmtData = pmtdata[:,[1,end]]
            if isempty(pmtData)
                continue
            end
            pdir = calculate_daily_pdir2(yyyymmdd, econ)
            if isempty(pdir)
                continue
            end
            out = value_Comparison(pmtData,pdir)
            if isempty(out)
                println("validation passed for econ $econ on date $yyyymmdd ")
            else
                result[econ][yyyymmdd] = out
#                failedEcon[iEcon,iday] = "$econ-$yyyymmdd"
                println("validation failed for econ $econ on date $yyyymmdd !!!")
            end
            if yyyymmdd == till_date
                break
            end
        end
        if isempty(result[econ])
            println("validation passed for econ $econ")
        end
    end
    return result
end


# This function is to validdate pdir for a single econ on a specified date since 20200101
function validate_daily_pdir2(yyyymmdd,econ)
    path = string(GC_VT["dataEnv"], "ProductionData\\Recent\\Daily\\")
    folder = readdir(path)
    date = folder[occursin.(string(yyyymmdd), folder)][1]
    pmtpath = string(path,date,"\\Products\\P17_Pdir2.0\\pdirNum_$econ.mat")
    pmtdata = vt_read_mat(pmtpath)
    pmtData = pmtdata[:,[1,end]]
    pdir = calculate_daily_pdir2(yyyymmdd, econ)
    out = value_Comparison(pmtData,pdir)
    if isempty(out)
        println("validation passed for econ $econ on date $yyyymmdd")
    else
        println("validation failed for econ $econ on date $yyyymmdd !!!")
    end
    return out
end
