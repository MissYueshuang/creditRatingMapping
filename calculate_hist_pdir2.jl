## Developed: 20200506

mpath = @__DIR__
cd(mpath)
k = findfirst("VT\\", mpath)
CommonToolPath = string(mpath[1:k[end]], "VT_common_tools\\Julia_common_tools\\")
include(string(CommonToolPath,"Include_All_Functions.jl"))
Include_All_Functions(CommonToolPath)
GC_VT = global_constants_VT()

function calculate_hist_pdir2(econ,yyyymm)
    histpath = string(GC_VT["dataEnv"], "ProductionData\\Historical\\$yyyymm\\Daily\\Products\\P2_Pd\\")
    resultDates = vt_read_mat(string(histpath,"$econ\\resultDates_$econ.mat"))
    resultFirms = vt_read_mat(string(histpath,"$econ\\resultFirms_$econ.mat"))
    resultPD = vt_read_mat(string(histpath,"$econ\\resultPD_$(econ)_12.mat"))
    boundry_path = raw"\\unicorn6\TeamData\VT\CRI_VT_MonitoringPDiR\\ratingInfo2.csv"
    boundry = CSV.read(boundry_path)
    if yyyymm/1000000 > 1 # add for revision date yyyymmdd
        resultDates = Int.(resultDates[:])
        if length(resultFirms) > 1
            resultFirms = Int.(resultFirms[:])
        else
            resultFirms = [Int.(resultFirms)]
        end
    end
    working_day = findWorkDay(resultDates[1],resultDates[end])
    lia, loc = ismember(resultDates,working_day)
    pdall = resultPD[lia,:]
    pdir = zeros(size(pdall))
    pdir .= NaN
    fisrt_valid_idx(x) = findall(.!isnan.(x))[1]
    fisrt_idx = mapslices(fisrt_valid_idx, pdall, dims=1)
    old_comps = []

    for iDate = minimum(fisrt_idx):length(working_day)
        io = working_day[iDate]
#        println("this is date $io")
        if iDate <= 10
            pd10day = pdall[1:iDate,:]
#            workingday = working_day[1:iDate]
#            workingday = sort!(workingday, rev = true)
        else
            pd10day = pdall[iDate-9:iDate,:]
#            workingday = working_day[iDate-9:iDate]
        end
        pdmean = mapslices(nanmean, pd10day, dims=1)
        pdir_daily = zeros(1,size(pdmean,2))
#        pdir_daily .= NaN
        init_Bound = boundry[2,:]
        upgrade_Bound = boundry[4,:]
        downgrade_Bound = boundry[5,:]
        for rBound = 1:length(init_Bound)
        #    pdir[pdmean.<=boundry[end+1-rBound]] = fill(length(boundry)+1-rBound,length(pdir[pdmean.<=boundry[end+1-rBound]]))
            pdir_daily[pdmean.<=float(init_Bound[end+1-rBound])] = fill(length(init_Bound)+1-rBound,length(pdir_daily[pdmean.<=init_Bound[end+1-rBound]]))
        end
        try
            pdir_daily[isnan.(pd10day[end,:])] .= NaN
        catch

        end

        if iDate > minimum(fisrt_idx)
            idx = findall(.!isnan.(pdir_daily))
            cols = [i[2] for i in idx]
            comps = resultFirms[cols]
            if isempty(comps)
                continue
            end
            lia, loc = ismember(comps,old_comps)
            cols = cols[lia]

            if !isempty(cols)
                last_pdir = pdir[1:(iDate-1),cols]
                find_last(x) = [findall(.!isnan.(x))[end], x[findall(.!isnan.(x))[end]]]
                last_pdir = mapslices(find_last, last_pdir, dims=1)
                last_pdir = transpose(last_pdir)
                temp_date = working_day[Int.(last_pdir[:,1])]
                last_pdir[:,1] = temp_date
                last_pdir = [comps[lia] last_pdir]

                pdir_daily_valid = transpose(pdir_daily[:,cols])
                pdmean_valid = copy(pdir_daily_valid)
                pdmean_valid[:] .= pdmean[cols]
                difference = last_pdir[:,end].-pdir_daily_valid

                pdir_unchanged = pdir_daily_valid[difference .== 0]
                idx_down = difference .< 0
                pdir_down = pdir_daily_valid[idx_down]
                for rBound = 1:length(downgrade_Bound)
                    pdir_down[pdmean_valid[idx_down].>=float(downgrade_Bound[rBound])] = fill(rBound,length(pdir_down[pdmean_valid[idx_down].>=downgrade_Bound[rBound]]))
                end
                idx_up = difference .> 0
                pdir_up = pdir_daily_valid[idx_up]
                for rBound = 1:length(upgrade_Bound)
                    pdir_up[pdmean_valid[idx_up].<=float(upgrade_Bound[end+1-rBound])] = fill(length(upgrade_Bound)-rBound+1,length(pdir_up[pdmean_valid[idx_up].<=upgrade_Bound[end+1-rBound]]))
                end
                pdir_daily_valid[difference .== 0] = pdir_unchanged
                pdir_daily_valid[idx_up] = pdir_up
                pdir_daily_valid[idx_down] = pdir_down

                pdir_daily[idx[lia]] = pdir_daily_valid
                pdir[iDate,:] = pdir_daily
                append!(old_comps,comps[:])
                old_comps = unique(old_comps)
            else
                append!(old_comps,comps[:])
                old_comps = unique(old_comps)
                pdir[iDate,:] = pdir_daily
            end
        else
            idx = findall(.!isnan.(pdir_daily))
            cols = [i[2] for i in idx]
            comps = resultFirms[cols]
            append!(old_comps,comps[:])
            pdir[iDate,:] = pdir_daily
        end
    end

    return pdir

end
