## Developed: 20200505
# This function id used to calculate incremental pdir2.0

using CSV
using DataFrames

mpath = @__DIR__
cd(mpath)
k = findfirst("VT\\", mpath)
CommonToolPath = string(mpath[1:k[end]], "VT_common_tools\\Julia_common_tools\\")
include(string(CommonToolPath,"findNextWorkday.jl"))
include(string(CommonToolPath,"findPrevday.jl"))
include(string(CommonToolPath,"vt_read_mat.jl"))
include(string(CommonToolPath,"ismember.jl"))
include(string(CommonToolPath,"getExcluded.jl"))
include(string(CommonToolPath,"nanmean.jl"))
include(string(CommonToolPath,"rowSetdiff.jl"))

#Include_All_Functions(CommonToolPath)
GC_VT = global_constants_VT()

function calculate_daily_pdir2(yyyymmdd, econ)
    HistDate = GC_VT["HistData_Date"]
    YearlyRefreshMonth = GC_VT["YearlyRefreshMonth"]
    firstDailyDate = findNextWorkday(HistDate)

    path = string(GC_VT["dataEnv"], "ProductionData\\Recent\\Daily\\")
    boundry_path = string(@__DIR__, "\\ratingInfo2.csv")
    boundry = CSV.read(boundry_path) #  for a file without column names (header), specify datarow=1

    workingday = findPrevday(yyyymmdd, 100) # randomly assignt a number 100 to make sure there are at least 10 valid days
    push!(workingday,yyyymmdd)
    workingday = sort!(workingday, rev = true)
    global valid_date = []
    global pd10day = []
    global compList_lastday = []
    folder = readdir(path) # all daily folders
    sort!(folder)
    flag = 0
    for iDate = 1:length(workingday)
#        global workingday
        date = workingday[iDate]
        folder_name = folder[occursin.(string(date), folder)]
        pdpath = string(path, folder_name[1], "\\Products\\P2_Pd\\")
        pdfile = string(pdpath, "pd60h_", string(econ), ".mat")
        pd = vt_read_mat(pdfile)
#        global pd
        if pd[1,2]*10000+ pd[1,3]*100+pd[1,4] != date
            if iDate == 1
#                println("Today is non-trading day, pdir2.0 copied last day data directly!!")
                flag = 1
                break
            else
                continue
            end
        end
        pd = pd[sortperm(pd[:, 1]), :]
    #    pd1 = sort(pd, dims=1, by = x->x[1])
        if iDate == 1 ## if no pd today, then no pdir today
            global compList_lastday = pd[:,1]
            global pd10day = pd[:,[1,16]]
        else
            compList = pd[:,1]
            idx = ismember(compList,compList_lastday)[1]
            pos = ismember(compList,compList_lastday)[2]
            pdlday = [NaN for i=1:size(compList_lastday,1), j=1]
        #    pdlday = zeros(size(compList_lastday,1),1)
            pdlday[pos[idx]] = pd[idx,16]
            pd10day = [pd10day pdlday]
        end
        push!(valid_date,date)
        if size(pd10day)[2] == 11
            workingday = valid_date
            break
        end
    end

    if flag == 1
        pd10day = []
    end
    if isempty(pd10day) # today is not trading days
        if yyyymmdd == firstDailyDate
            pdir = lastPdir(econ, folder, path, workingday,YearlyRefreshMonth,firstDailyDate)
            date = workingday[1]
            folder_name = folder[occursin.(string(date), folder)]
            pdpath = string(path, folder_name[1], "\\Products\\P2_Pd\\")
            pdfile = string(pdpath, "pd60h_", string(econ), ".mat")
            pd = vt_read_mat(pdfile)
            loc, lia = ismember(pdir[:,1],pd[:,1].*1000)
            pdir = pdir[loc,:]
            pdir[:,1] = floor.(pdir[:,1]./1000)
            pdir = pdir[:,1:end .!=2]
        else
            folder_name = folder[occursin.(string(workingday[2]), folder)][1]
            pdirpath = string(path,folder_name,"\\Products\\P17_Pdir2.0\\pdirNum_",string(econ), ".mat")
            pdir = vt_read_mat(pdirpath)
            pdir = pdir[:,[1,end]]
        end
        return pdir
#        return pd10day
    else ### trading day, compute pdir
        ## compute initial pdir
        # consider ExcludedFirms
        ExcludeFirms = getExcluded(workingday[10],yyyymmdd)
        idx2 = ismember(ExcludeFirms[:,1],pd10day[:,1])[1]
        excluded = ExcludeFirms[idx2,:]
        excludeRow = []
        for iEx = 1:size(excluded)[1]
            if (workingday[1] < excluded[iEx,4]) | isnan(excluded[iEx,4]) # still excluded
                pd10day[pd10day[:,1].==excluded[iEx,1],2:end] .= NaN
                push!(excludeRow,findall(pd10day[:,1].==excluded[iEx,1])[1])
            else
                startday = findall(workingday.>=excluded[iEx,3])[end]
                stopday = findall(workingday.>=excluded[iEx,4])[end]
                pd10day[pd10day[:,1].==excluded[iEx,1],stopday+2:startday+1] .= NaN
            end
        end
    #    pd10day[isnan.(pd10day)] .= 0

    #    nanzeromean(x) = mean(filter( !iszero , x ))
        pdmean = mapslices(nanmean, pd10day[:,2:end], dims=2)
        idx = findall(isnan.(pdmean))
        pdmean[idx] .= 0
        pdir = zeros(size(pdmean,1),1)
        init_Bound = boundry[2,:]
        upgrade_Bound = boundry[4,:]
        downgrade_Bound = boundry[5,:]
        for rBound = 1:length(init_Bound)
            pdir[pdmean.<=float(init_Bound[end+1-rBound])] = fill(length(init_Bound)+1-rBound,length(pdir[pdmean.<=init_Bound[end+1-rBound]]))
        end
        try
            pdir[isnan.(pd10day[:,end])] = NaN
        catch

        end
        # apply buffer zone, compare with last valid pdir value
        ##### here is load directly from processing
        # folder_name = folder[occursin.(string(yyyymmdd), readdir(path))][1]
        # last_pdir_file = string(path,folder_name,"\\Processing\\P17_Pdir2.0\\lastRating_",string(econ), ".mat")
        # last_pdir = vt_read_mat(last_pdir_file)
        last_pdir = lastPdir(econ, folder, path, workingday,YearlyRefreshMonth,firstDailyDate)
        lia,loc = ismember(last_pdir[:,1],pd10day[:,1].*1000)
        last_valid_pdir = last_pdir[lia,:]
        last_valid_pdir = last_valid_pdir[sortperm(last_valid_pdir[:, 1]), :]

        dif, diffrow = rowSetdiff(pd10day[:,1].*1000,last_valid_pdir[:,1]) # if exist new company
        pdir_new_comp = pdir[diffrow]
        if !isempty(diffrow)
            pdir = pdir[setdiff(1:end,diffrow),:] # ignore new company first
            pdmean = pdmean[setdiff(1:end,diffrow),:]
        end
        difference = last_valid_pdir[:,end].-pdir
        pdir_unchanged = pdir[difference .== 0]
        idx_down = difference .< 0
        pdir_down = pdir[idx_down]
        for rBound = 1:length(downgrade_Bound)
            pdir_down[pdmean[idx_down].>=float(downgrade_Bound[rBound])] = fill(rBound,length(pdir_down[pdmean[idx_down].>=downgrade_Bound[rBound]]))
        end
        idx_up = difference .> 0
        pdir_up = pdir[idx_up]
        for rBound = 1:length(upgrade_Bound)
            pdir_up[pdmean[idx_up].<=float(upgrade_Bound[end+1-rBound])] = fill(length(upgrade_Bound)-rBound+1,length(pdir_up[pdmean[idx_up].<=upgrade_Bound[end+1-rBound]]))
        end
        pdir[difference .== 0] = pdir_unchanged
        pdir[idx_up] = pdir_up
        pdir[idx_down] = pdir_down

        pdir = [pd10day[setdiff(1:end,diffrow),1] pdir]
        new_pdir = [pd10day[diffrow,1] pdir_new_comp]
        if !isempty(new_pdir)
            pdir = [pdir;new_pdir]
        end
        pdir = pdir[sortperm(pdir[:, 1]), :]
        pdir[excludeRow,2] .= 0
    end
    return pdir

end

function lastPdir(econ, folder, path, workingday, YearlyRefreshMonth, firstDailyDate)
    hist_path = string(GC_VT["dataEnv"],"ProductionData\\Historical\\$YearlyRefreshMonth\\Daily\\Products\\P17_Pdir2.0\\")
    hist_pdir = vt_read_mat(string(hist_path,"resultRating_$econ.mat"))
    hist_comp = vt_read_mat(string(hist_path,"resultFirms_$econ.mat"))
    hist_date = vt_read_mat(string(hist_path,"resultDates_$econ.mat"))

    find_last(x) = [findall(.!isnan.(x))[end], x[findall(.!isnan.(x))[end]]]
    hist_pdir = mapslices(find_last, hist_pdir, dims=1)
    hist_pdir = transpose(hist_pdir)
    hist_date = hist_date[Int.(hist_pdir[:,1])]
    hist_pdir[:,1] = hist_date
    hist_pdir = [hist_comp hist_pdir] # [comp, date, pdir]

    baseyear = Int(floor(firstDailyDate/10000))
    folder_new = folder[(occursin.(string(baseyear),folder)) .| (occursin.(string(baseyear+1),folder))]
#    folder_new = folder[occursin.(string(Int(floor(firstDailyDate/10000))),folder)]
    sort!(folder_new)
    idx = findall(startswith.(folder_new,string(firstDailyDate)))[1]
    folder_new = folder_new[idx:end]

    lastday = workingday[2]
    if floor(lastday[1]/100) == YearlyRefreshMonth
        return hist_pdir
    end
    for date = folder_new
        if !isfile(string(path,date,"\\Products\\P17_Pdir2.0\\pdirNum_$econ.mat"))
            if parse(Int64,date[1:8]) == lastday
                return hist_pdir
            else
                continue
            end
        else
            daily_pdir = vt_read_mat(string(path,date,"\\Products\\P17_Pdir2.0\\pdirNum_$econ.mat"))
        end
        idx = findall(daily_pdir[:,5].==0)
        daily_pdir = daily_pdir[setdiff(1:end,idx),:]
        lia,loc = ismember(hist_pdir[:,1],daily_pdir[:,1].*1000)
        hist_pdir[lia,3] = daily_pdir[loc[loc.!=0],5]
        hist_pdir[lia,2] = daily_pdir[loc[loc.!=0],2]*10000+daily_pdir[loc[loc.!=0],3]*100+daily_pdir[loc[loc.!=0],4]
        dif,diffrow = rowSetdiff(daily_pdir[:,1]*1000,hist_pdir[:,1])
        if !isempty(dif)
            firstDate = daily_pdir[diffrow,2]*10000+daily_pdir[diffrow,3]*100+daily_pdir[diffrow,4]
            first_pdir = daily_pdir[diffrow,5]
            new_comp = daily_pdir[diffrow,1].*1000
            new_data = [new_comp firstDate first_pdir]
            hist_pdir = [hist_pdir;new_data]
        end
        if parse(Int64,date[1:8]) == lastday
            hist_pdir[:,1] .= Int.(hist_pdir[:,1])
            return hist_pdir
        end
    end
end
