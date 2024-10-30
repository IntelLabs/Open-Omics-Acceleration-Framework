if [ -z "$1" ]; then
    echo "Provide the protein name"
    exit 1
fi
url="https://zenodo.org/records/4031961/files/data.zip?download=1"
download_dir="./data_new"
target_folder="$1"  
mkdir -p "$download_dir"
echo "Downloading data.zip..."
wget -O "$download_dir/data.zip" "$url"
echo "Unzipping data.zip..."
unzip "$download_dir/data.zip" -d "$download_dir"
if [ -d "$download_dir/data" ]; then
    echo "Nested 'data' folder found."
    if [ -d "$download_dir/data/$target_folder" ]; then
        mv "$download_dir/data/$target_folder" ./
        echo "$target_folder folder successfully moved to the current directory."
    else
        echo "$target_folder folder not found inside nested 'data' directory."
    fi
    rm -rf "$download_dir/data"
fi
rm -rf "$download_dir"
echo "Download and extraction complete. '$target_folder' folder is now in the current directory."
