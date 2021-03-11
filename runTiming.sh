runTiming () {
    
    if [ "$(uname)" == "Darwin" ]; then
        sed -e "s/^TYPE[[:blank:]].*$/TYPE = ${1}/" -i '' Makefile.in
        sed -e "s/^SUBTYPE[[:blank:]].*$/SUBTYPE = ${2}/" -i '' Makefile.in
    else
        sed -i "s/^TYPE\s.*$/TYPE = ${1}/" Makefile.in
        sed -i "s/^SUBTYPE\s.*$/SUBTYPE = ${2}/" Makefile.in
    fi

    echo "Compiling ${1} ${2}"
    make -s clean
    make -s
    echo "Running ${1} ${2}"
    ./disco
    mv times.log times.${1}.${2}.log
    mv grid.txt grid.${1}.${2}.txt
}

runTiming 0 0
runTiming 1 0
runTiming 1 1
