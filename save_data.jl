using JSON

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