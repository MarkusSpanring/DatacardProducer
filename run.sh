make -B

if [ ! -d "histos" ]; then
    mkdir -p histos;
fi
echo $1
./makeHTTDatacards $1
