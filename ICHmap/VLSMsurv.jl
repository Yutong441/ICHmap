using Survival, StatsModels, CSV, DataFrames, GLM, NIfTI, ArgParse
using MultipleTesting, HypothesisTests


function run_cox(df, xvar, outcome, covars)
    # `df`: dataframe, must contain outcome"_time" and outcome"_bool" and the
    # covariates
    # `xvar`: a vector, same length as the number of rows in `df`
    df.x = xvar
    no_na = dropmissing(df, [outcome*"_time", outcome*"_bool"])
    no_na.y = EventTime.(no_na[:, outcome*"_time"], no_na[:, outcome*"_bool"] .== 1)

    all_vars = [Term(Symbol("x"))]
    for i in 1:length(covars)
        push!(all_vars, Term(Symbol(covars[i])))
    end
    model = coxph(FormulaTerm(Term(Symbol("y")), Tuple(all_vars)), no_na)

    ctab = coeftable(model)
    ctab = DataFrame(ctab)

    # create a vector of voxels and covariates
    M = length(covars)
    var_stat = zeros(2 + 2*M)
    var_stat[1] = ctab[(ctab.Name .== "x"), "Pr(>|z|)"][1]
    var_stat[2] = exp(ctab[(ctab.Name .== "x"), "Estimate"][1])
    for i in 1:M
        var_stat[1+i*2] = ctab[(ctab.Name .== covars[i]), "Pr(>|z|)"][1]
        var_stat[2+i*2] = exp(ctab[(ctab.Name .== covars[i]), "Estimate"][1])
    end
    return var_stat
end


function sum_map(all_files, threshold, save_path) 
    # find out which voxels happen at a frequency above the threshold
    one_img = niread(all_files[1])
    sum_img = one_img.raw
    for i in 2:length(all_files)
        img = niread(all_files[i])
        sum_img = sum_img + img.raw
    end
    save_nii(sum_img, one_img, save_path)
    proceed_map = sum_img .> threshold*length(all_files)
    return proceed_map
end


function obtain_coord(proceed_map)
    dimen = size(proceed_map)
    N = floor(Int, sum(proceed_map))
    coord = zeros(Int, (N, 3))
    index = 1
    for k in 1:dimen[3]
        for j in 1:dimen[2]
            for i in 1:dimen[1]
                if (proceed_map[i, j, k] == true)
                    coord[index, 1] = i
                    coord[index, 2] = j
                    coord[index, 3] = k
                    index = index + 1
                end
            end
        end
    end
    return coord
end

function obtain_matrix(all_files, coord)
    N = length(all_files)
    L = size(coord)[1]
    all_mat = zeros(N, L)
    for i in 1:N
        img = niread(all_files[i])
        for j in 1:L
            x = coord[j, 1]
            y = coord[j, 2]
            z = coord[j, 3]
            all_mat[i, j] = img.raw[x, y, z]
        end
    end
    return all_mat
end

function run_cox_all(all_df, all_mat, outcome, covar)
    N = size(all_mat)[1]
    L = size(all_mat)[2]
    M = length(covar)
    out_stat = zeros(L, 2+2*M)
    for i in 1:L
        out = run_cox(all_df, all_mat[:, i], outcome, covar)
        out_stat[i, :] = out
    end

    for i in 1:(M+1)
        out_stat[:, 2*i-1] = adjust(out_stat[:, 2*i-1], BenjaminiHochberg())
    end
    return out_stat
end


function coord_to_map(vec, coord, img_size)
    # convert tables of p values or other statistical outcomes into maps
    N = size(vec)[1]
    TY = size(vec)[2]
    img = ones(img_size[1], img_size[2], img_size[3], TY)
    for j in 1:TY
        for i in 1:N
            x = coord[i, 1]
            y = coord[i, 2]
            z = coord[i, 3]
            img[x, y, z, j] = vec[i, j]
        end
    end
    return img
end

function save_nii(arr, nifti_ob, save_path)
    arr_nib = NIVolume(nifti_ob.header, arr)
    niwrite(save_path, arr_nib)
end

function comp_time(all_df, all_mat, pval_mat, outcome)
    # all_df: demographics dataframe
    # pval_mat: p value array from running `run_cox_all`
    # all_mat: binary array containing lesion status in every voxel in every
    # patient
    N = size(all_mat)[1]
    L = size(all_mat)[2]
    duration = all_df[:, outcome*"_time"]  # N
    status = all_df[:, outcome*"_bool"]

    time_mat = zeros(L, 2)
    for i in 1:L
        # only test for the significant voxels
        if pval_mat[i, 1] < 0.1
            lesion_status = all_mat[:, i]
            # only include patients who developed a particular outcome
            lesion_status = lesion_status[status .== 1]
            duration_vox = duration[status .== 1]
            lesion_time = duration_vox[lesion_status .== 1]
            nolesion_time = duration_vox[lesion_status .== 0]
            mod = MannWhitneyUTest(lesion_time, nolesion_time)
            time_mat[i, 1] = pvalue(mod)
            time_mat[i, 2] = sum(lesion_time)/length(lesion_time)/sum(nolesion_time)*length(nolesion_time)
        else
            time_mat[i, 1] = 1
            time_mat[i, 2] = 0
        end
    end
    time_mat[:, 1] = adjust(time_mat[:, 1], BenjaminiHochberg())
    return time_mat
end

function vox_wise_cox(root, outcome, save_dir, covar, vox_thres)
    # save_dir: save 2 images: stat_img.nii.gz (map of hazard ratio),
    # p_val.nii.gz (p value map after FDR correction)
    println("reading csv")
    df_path = root*"/"*outcome*".csv"
    all_df = CSV.read(df_path, DataFrame)
    all_paths = all_df.paths

    println("loading maps")
    proceed_map = sum_map(all_paths, vox_thres, save_dir*"/freq_map.nii.gz")
    coord_one = obtain_coord(proceed_map);
    all_mat = obtain_matrix(all_paths, coord_one);

    println("voxel wise Cox regression")
    out_mat = run_cox_all(all_df, all_mat, outcome, covar)

    println("saving output")
    one_img = niread(all_paths[1])
    out_maps = coord_to_map(out_mat, coord_one, size(one_img))

    M = length(covar)
    save_nii(out_maps[:, :, :, 1], one_img, save_dir*"/pval_img.nii.gz")
    save_nii(out_maps[:, :, :, 2], one_img, save_dir*"/stat_img.nii.gz")
    # for i in 1:M
    #     save_nii(out_maps[:, :, :, 2*i+1], one_img,
    #              save_dir*"/"*covar[i]*"_pval.nii.gz")
    #     save_nii(out_maps[:, :, :, 2*i+2], one_img,
    #              save_dir*"/"*covar[i]*"_stat.nii.gz")
    # end

    println("timing analysis")
    time_mat = comp_time(all_df, all_mat, out_mat, outcome)
    out_time = coord_to_map(time_mat, coord_one, size(one_img))
    save_nii(out_time[:, :, :, 1], one_img, save_dir*"/pval_time.nii.gz")
    save_nii(out_time[:, :, :, 2], one_img, save_dir*"/stat_time.nii.gz")
end


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--root"
            help = "root to the image and behavior dataframes"
            required = true
        "--outcome"
            help = "which outcome to predict"
            required = true
        "--save_dir"
            help = "where to save the data"
            required = true
        "--covar"
            help = "covariates (comma separated)"
            default = "Age,Sex,ICH_vol"
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    if !isdir(parsed_args["save_dir"])
        mkdir(parsed_args["save_dir"])
    end
    covar = collect(eachsplit(parsed_args["covar"], ","))
    vox_wise_cox(parsed_args["root"], parsed_args["outcome"],
                 parsed_args["save_dir"], covar, 0.025)
end

@time main()
