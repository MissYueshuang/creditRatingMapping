## Developed: 20200506
using CSV

mpath = @__DIR__
cd(mpath)
k = findfirst("VT\\", mpath)
CommonToolPath = string(mpath[1:k[end]], "VT_common_tools\\Julia_common_tools\\")
include(string(CommonToolPath,"findNextWorkday.jl"))
include(string(CommonToolPath,"value_Comparison.jl"))
include(string(CommonToolPath,"global_constants_VT.jl"))
include(string(CommonToolPath,"getAllEcon.jl"))
# Include_All_Functions(CommonToolPath)
GC_VT = global_constants_VT()
include(string(GC_VT["vtEnv"],"CRI_VT_MonitoringPDiR\\PDiR2.0\\calculate_daily_pdir2.jl"))

## This function is to validate incremental pdir2.0 on a specific daily date
function validate_daily_PDiR2(numPDdate, Econ=[])
    HistDate = GC_VT["HistData_Date"]
    YearlyRefreshMonth = GC_VT["YearlyRefreshMonth"]
    firstDailyDate = findNextWorkday(HistDate)

    path = string(GC_VT["dataEnv"], "ProductionData\\Recent\\Daily\\")
    folder = readdir(path)
    if isempty(Econ)
        econList = getAllEcons()
    else
        econList = Econ
    end
    result = Dict()
    # n = length(econList)
    # folder_new = folder[occursin.("2020",folder)]
    baseyear = Int(floor(firstDailyDate/10000))
    folder_new = folder[(occursin.(string(baseyear),folder)) .| (occursin.(string(baseyear+1),folder))]
    date = folder_new[startswith.(folder_new,string(numPDdate))][1]
    # till_folder_name = folder[occursin.(string(till_date), folder)][1]
    # m = findall(folder_new.==till_folder_name)[1]
    # failedEcon = SharedArray{Float64, 2}((n,m))

    for econ = econList
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
            result["Econ$econ"] = out
#                failedEcon[iEcon,iday] = "$econ-$yyyymmdd"
            println("validation failed for econ $econ on date $yyyymmdd !!!")
        end
    end
    return result
end
