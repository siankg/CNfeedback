#!/bin/sh

# Set the following directories:
dir_raw= #location of raw files
dir_processed= #location of processed files

# Merge time:
experiments=("1pctCO2 1pctCO2-bgc 1pctCO2-rad 1pctCO2Ndep 1pctCO2Ndep-bgc")
for experiment in $experiments
do
    models=("ACCESS-ESM1-5 BCC-CSM2-MR CanESM5 CESM2 CNRM-ESM2-1 GFDL-ESM4 IPSL-CM6A-LR MIROC-ES2L MPI-ESM1-2-LR NorESM2-LM UKESM1-0-LL")
    for model in $models
    do

        # netAtmosLandCO2Flux
        if [ "${model}" = "BCC-CSM2-MR" ]; then
            cdo -s -L -yearmean -selyear,1850/1989 -mergetime ${dir_raw}/nep*${model}*${experiment}*.nc ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
            ncrename -v nep,netAtmosLandCO2Flux ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
        elif [ "${model}" = "ACCESS-ESM1-5" ]  || [ "${model}" = "GFDL-ESM4" ] || [ "${model}" = "NorESM2-LM" ]; then
            cdo -s -L -yearmean -seltimestep,1/1680 -mergetime ${dir_raw}/nbp*${model}*${experiment}*.nc ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
            ncrename -v nbp,netAtmosLandCO2Flux ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
        elif [ "${model}" = "IPSL-CM6A-LR" ]; then
            cdo -s -L -yearmean -selyear,1850/1989 -mergetime ${dir_raw}/nbp*${model}*${experiment}*.nc ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
            ncrename -v nbp,netAtmosLandCO2Flux ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
        elif [ "${model}" = "CESM2" ]; then
            cdo -s -L -yearmean -seltimestep,1/1680 -mergetime ${dir_raw}/netAtmosLandCO2Flux*${model}*${experiment}*.nc ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
        else
            cdo -s -L -yearmean -selyear,1850/1989 -mergetime ${dir_raw}/netAtmosLandCO2Flux*${model}*${experiment}*.nc ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
        fi

        # fgco2 and tas
        if [ "${model}" = "ACCESS-ESM1-5" ]  || [ "${model}" = "CESM2" ] || [ "${model}" = "GFDL-ESM4" ] || [ "${model}" = "NorESM2-LM" ]; then
            cdo -s -L -yearmean -seltimestep,1/1680 -mergetime ${dir_raw}/fgco2*${model}*${experiment}*.nc ${dir_processed}/fgco2_${model}_${experiment}.nc
            cdo -s -L -yearmean -seltimestep,1/1680 -mergetime ${dir_raw}/tas*${model}*${experiment}*.nc ${dir_processed}/tas_${model}_${experiment}.nc
        else
            cdo -s -L -yearmean -selyear,1850/1989 -mergetime ${dir_raw}/fgco2*${model}*${experiment}*.nc ${dir_processed}/fgco2_${model}_${experiment}.nc
            cdo -s -L -yearmean -selyear,1850/1989 -mergetime ${dir_raw}/tas*${model}*${experiment}*.nc ${dir_processed}/tas_${model}_${experiment}.nc
        fi

        variables=("fNnetmin fBNF cVeg nVeg fN2O")
        for variable in $variables
        do
            cdo -s -L -yearmean -selyear,1850/1989 -mergetime ${dir_raw}/${variable}*${model}*${experiment}*.nc ${dir_processed}/${variable}_${model}_${experiment}.nc
        done

        variables=("areacella areacello sftlf sftof")
        for variable in $variables
        do
            cp ${dir_raw}/${variable}*${model}*.nc ${variable}_${model}.nc
        done

    done
done

# Calculate global totals, means, and zonal sums for plotting:
experiments=("1pctCO2 1pctCO2-bgc 1pctCO2-rad 1pctCO2Ndep 1pctCO2Ndep-bgc")
for experiment in $experiments
do
    models=("ACCESS-ESM1-5 BCC-CSM2-MR CanESM5 CESM2 CNRM-ESM2-1 GFDL-ESM4 IPSL-CM6A-LR MIROC-ES2L MPI-ESM1-2-LR NorESM2-LM UKESM1-0-LL")
    for model in $models
    do

                file=${dir_processed}/fgco2_${model}_${experiment}.nc
                areafile=${dir_processed}/areacello_${model}.nc
                fractionfile=${dir_processed}/sftof_${model}.nc
                cdo -s -L -fldsum -mul -mul $file $areafile -divc,100 $fractionfile ${dir_processed}/fgco2_${model}_${experiment}_total.nc
                cdo -s -remapbil,r360x180  -mul -mul $file $areafile -divc,100 $fractionfile temp.nc
                cdo -s -timsum -zonsum temp.nc  ${dir_processed}/fgco2_${model}_${experiment}_zonsum_timsum.nc
                rm temp.nc

                file=${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}.nc
                areafile=${dir_processed}/areacella_${model}.nc
                fractionfile=${dir_processed}/sftlf_${model}.nc
                cdo -s -L -fldsum -mul -mul $file $areafile -divc,100 $fractionfile ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}_total.nc
                cdo -s -remapbil,r360x180  -mul -mul $file $areafile -divc,100 $fractionfile temp.nc
                cdo -s -timsum -zonsum temp.nc  ${dir_processed}/netAtmosLandCO2Flux_${model}_${experiment}_zonsum_timsum.nc

                file=${dir_processed}/tas_${model}_${experiment}.nc
                areafile=${dir_processed}/areacella_${model}.nc
                cdo -s -L -div -fldsum -mul $file $areafile -fldsum $areafile ${dir_processed}/tas_${model}_${experiment}_mean.nc
                cdo -s -enssum -seltimestep,140 -zonmean $file -mulc,-1 -seltimestep,1 -zonmean $file ${dir_processed}/tas_${model}_${experiment}_zonmean_timdiff.nc

                variables=("fNnetmin fBNF cVeg nVeg fN2O")
                for variable in $variables
                do
                    file=${dir_processed}/${variable}_${model}_${experiment}.nc
                    areafile=${dir_processed}/areacella_${model}.nc
                    fractionfile=${dir_processed}/sftlf_${model}.nc
                    cdo -s -L -fldsum -mul -mul $file $areafile -divc,100 $fractionfile ${dir_processed}/${variable}_${model}_${experiment}_total.nc
                done

    done
done

# Calculate total and zonal sums of N deposition for plotting:
cdo -s -yearmean -enssum drynhx_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc \
    drynoy_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc \
    wetnhx_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc \
    wetnoy_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc \
    Ndep_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc

cdo -s -L -selyear,1850/1989 -fldsum -mul Ndep_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc \
    -gridarea Ndep_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc \
    Ndep_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012_total.nc

cdo -s -L -mul Ndep_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc -gridarea Ndep_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc temp.nc
cdo -s -enssum -zonsum -selyear,1989 temp.nc -mulc,-1 -zonsum -selyear,1850 temp.nc Ndep_zonsum_timdiff.nc
cdo -s -zonsum -timsum -selyear,1850/1989 temp.nc Ndep_zonsum_timsum.nc
rm temp.nc