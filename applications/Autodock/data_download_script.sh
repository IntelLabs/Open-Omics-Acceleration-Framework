repo_url="https://github.com/emascarenhas/AutoDock-GPU.git"
data_dir="1ac8"
mkdir -p "$data_dir"
cd "$data_dir"
git init
git remote add origin "$repo_url"
git config core.sparseCheckout true
echo "input/1ac8/derived/*" >> .git/info/sparse-checkout
git pull origin develop
mv input/1ac8/derived/* .
rm -rf .git input
