url="https://zenodo.org/records/4031961/files/data.zip?download=1"
download_dir="./data_original"
target_folder="$1"  
mkdir -p "$download_dir"
echo "Downloading data.zip..."
wget -O "$download_dir/data.zip" "$url"
echo "Unzipping data.zip..."
unzip "$download_dir/data.zip" -d "$download_dir"
if [ -d "$download_dir/data" ]; then
    echo "Nested 'data' folder found."
    if [ -d "$download_dir/data/$target_folder" ]; then
        cp -r "$download_dir/data/$target_folder" ./
        echo "$target_folder folder successfully copied to the current directory."
    else
        echo "$target_folder folder not found inside nested 'data' directory."
    fi
else
  echo "'data' folder not found in the extracted files."
fi
rm -f "$download_dir/data.zip"
echo "'$target_folder' folder is now copied to the current directory, and the complete original data is retained in '$download_dir/data'."
