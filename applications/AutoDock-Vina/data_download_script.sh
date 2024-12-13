url="https://zenodo.org/records/4031961/files/data.zip?download=1"
download_dir="./data_original"
target_folder="$1"
if [ ! -d "$download_dir/data" ]; then
    echo "Downloading data.zip..."
    mkdir -p "$download_dir"
    wget -O "$download_dir/data.zip" "$url"

    echo "Unzipping data.zip..."
    unzip "$download_dir/data.zip" -d "$download_dir"
    rm -f "$download_dir/data.zip"

    echo "Data downloaded and extracted to $download_dir/data"
else
    echo "Data already exists in $download_dir/data. Skipping download and extraction."
fi
if [ -d "$target_folder" ]; then
    echo "The folder '$target_folder' already exists in the current directory. Skipping copy."
else
    if [ -d "$download_dir/data/$target_folder" ]; then
        cp -r "$download_dir/data/$target_folder" ./
        echo "$target_folder folder successfully copied to the current directory."
    else
        echo "$target_folder folder not found inside '$download_dir/data'."
    fi
fi
echo "'$target_folder' folder is now available in the current directory."
