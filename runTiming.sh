runTiming () {
    
    if [ "$(uname)" == "Darwin" ]; then
        sed -e "s/^TYPE[[:blank:]].*$/TYPE = ${1}/" -i '' Makefile.in
        sed -e "s/^SUBTYPE[[:blank:]].*$/SUBTYPE = ${2}/" -i '' Makefile.in
    else
        sed -i "s/^TYPE\s.*$/TYPE = ${1}/" Makefile.in
        sed -i "s/^SUBTYPE\s.*$/SUBTYPE = ${2}/" Makefile.in
    fi

    make clean
    make
    ./disco
    mv times.log times.${1}.${2}.log
}

runTiming 0 0
runTiming 1 0
runTiming 1 1
