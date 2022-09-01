runAll () {
    
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
    ./disco > out.${1}${2}.txt
}

runAll 0 0
runAll 1 0
runAll 1 1
runAll 2 0
