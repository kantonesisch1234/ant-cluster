using JLD2, FileIO, JSON

function save_data_to_json(data,filename)
    data_str = JSON.json(data)
    open(filename,"w") do f
        write(f,data_str)
    end
end

function read_data_from_json(filename)
    data = JSON.parsefile(filename)
    return data
end

function save_data_to_jld2(data,filename)
	save(filename, Dict("data" => data))
end

function read_data_from_jld2(filename)
	data = load(filename, "data")
	return data
end

function get_file_extension(filename)
    try
        dot_idx = findlast(isequal('.'),filename)
    catch
        return 0
    end
    return filename[dot_idx:end]
end

function get_filename_wo_extension(filename)
    dot_idx = findlast(isequal('.'),filename)
    return filename[1:dot_idx-1]
end