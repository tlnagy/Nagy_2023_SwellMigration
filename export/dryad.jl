using HTTP
using JSON
using URIs

######################################
# Get access token
######################################

root = "https://datadryad.org"
doi = URIs.escapeuri("doi:10.7272/Q6NS0S5N")

app_id = ENV["AppID"]
secret = ENV["AppSecret"]

payload = Dict(
    "grant_type" => "client_credentials",
    "client_id" => app_id,
    "client_secret" => secret
)

response = HTTP.request("POST", "$(root)/oauth/token", ["Content-Type" => "application/json"], JSON.json(payload))
response = JSON.parse(String(response.body))
access_token = response["access_token"]

headers = Dict(
    "Accept" => "application/json",
    "Content-Type" => "application/json",
    "Authorization" => "Bearer $access_token"
)

######################################
# Get versions
######################################

resp = HTTP.get("$(root)/api/v2/datasets/$(doi)/versions", headers)
s = JSON.parse(String(resp.body))

version = s["_embedded"]["stash:versions"][end]
versionpath = version["_links"]["stash:files"]["href"]
@info "Downloading from version $(version["versionNumber"]), released $(version["lastModificationDate"])"

"""
Searches through DataDryads file index for a given filename, returns either a
dictionary or nothing if the file isn't found
"""
function getfilepath(root, versionpath, filename, headers)
    i = 1
    fileinfo = nothing
    while true
        files = HTTP.get(root * versionpath * "?page=$i", headers)
        filedict = JSON.parse(String(files.body))

        filelist = filedict["_embedded"]["stash:files"]
        foundfile = findfirst(x->occursin(filename, x["path"]), filelist)
        if foundfile !== nothing
            fileinfo = filelist[foundfile]
            break
        end

        if "next" ∉ keys(filedict["_links"])
            @info "File not found"
            return nothing
        end
        i += 1
    end
    return fileinfo
end

function download(root, versionpath, dst, filename, headers)
    fileinfo = getfilepath(root, versionpath, filename, headers)
    fileurl = fileinfo["_links"]["stash:file-download"]["href"]
    mkpath(dst)
    HTTP.download(root * fileurl, joinpath(dst, filename), headers; update_period = 10)
end


# download(root, versionpath, "data", "fxm_uncaging_augmented.csv", headers)


