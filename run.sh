make -B

if [ ! -d "histos" ]; then
    mkdir -p histos;
fi
echo $1
./runFile $1
