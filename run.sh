make -B

if [ ! -d "histos" ]; then
    mkdir -p histos;
fi

./runFile test
