#!/bin/bash

testmode=$(echo "$2" | tr '[:upper:]' '[:lower:]')

if [ "$testmode" = "archer" ] ; then
    declare -a compiler_targets=( "GNU" "Intel" "Cray" )
elif [ "$testmode" = "gnu" ] ; then
    testmode="local"
    declare -a compiler_targets=( "GNU" )
elif [ "$testmode" = "cray" ] ; then
    testmode="local"
    declare -a compiler_targets=( "Cray" )
elif [ "$testmode" = "intel" ] ; then
    testmode="local"
    declare -a compiler_targets=( "Intel" )
else
    testmode="local"
    declare -a compiler_targets=( "GNU" )
fi

declare -a test_names
declare -a test_monccomponents
declare -a test_unittests
declare -a test_makefiledir
declare -a test_testurl

unitTests=0
compilerPass=1
entirePass=1

resultsFile="unitruns.txt"
relativeFruitLocation="misc/fruit"

main() {
    parseTests $1

    if [ "$testmode" = "archer" ] ; then
        svn update
        module unload PrgEnv-cray
        module unload PrgEnv-gnu
        module unload PrgEnv-intel
    fi

    echo "[Start]$(date +%d/%m/%Y--%T)" > $resultsFile
    echo "Compiler targets:"${compiler_targets[@]} >> $resultsFile
    echo "Test Categories:"${unitTests} >> $resultsFile
    echo "[SVN]:"`svnversion | sed -e 's/M//g'` >> $resultsFile
    for i in "${compiler_targets[@]}"
    do
        runCompilerTests $i
    done
    if [ "$entirePass" = "1" ] ; then
        echo "[Entire]Pass" >> $resultsFile
        if [ "$testmode" = "local" ] ; then echo "[Entire]Pass" ; fi
    else
        echo "[Entire]Fail" >> $resultsFile
        if [ "$testmode" = "local" ] ; then echo "[Entire]Fail" ; fi
    fi
    echo "[End]$(date +%d/%m/%Y--%T)" >> $resultsFile
}

runCompilerTests() {
    compilerPass=1
    if [ "$testmode" = "archer" ] ; then
        eval "module load PrgEnv-"$(echo "$1" | tr '[:upper:]' '[:lower:]')
    fi
    compileFruit $1
    runAllTests $1
    if [ "$testmode" = "archer" ] ; then
        module load netcdf
        eval "make $1 EXEC_NAME=test-build-$1 &> tempcompile"
        if [ -f test-build-$1 ]; then
            echo "[Entire Build]"$1":Pass" >> $resultsFile
        else
            echo "[Entire Build]"$1":Fail" >> $resultsFile
            compilerPass=0
            entirePass=0
        fi
        cat tempcompile >> $resultsFile
        echo "--End Message--" >> $resultsFile
        rm tempcompile
        eval "rm -f test-build-$1"
        module unload netcdf
        eval "module unload PrgEnv-"$(echo "$1" | tr '[:upper:]' '[:lower:]')
    fi
    if [ "$compilerPass" = "1" ] ; then
        echo "[Compiler]"$1":Pass" >> $resultsFile
        if [ "$testmode" = "local" ] ; then echo "[Compiler]"$1":Pass" ; fi
    else
        echo "[Compiler]"$1":Fail" >> $resultsFile
        if [ "$testmode" = "local" ] ; then echo "[Compiler]"$1":Fail" ; fi
    fi
}

compileFruit() {
    pushd $relativeFruitLocation > /dev/null
    make $1 > /dev/null
    popd > /dev/null
}

runAllTests() {
    for (( i=0; i<$unitTests; i++ ))
    do
        runspecifictest "${test_monccomponents[i]}" "${test_unittests[i]}" ${test_names[i]} $1 "${test_makefiledir[i]}" "${test_testurl[i]}"
    done
}

parseTests() {
    file=$1

    while IFS=: read col1 col2 col3 col4 col5
    do
        test_names[unitTests]="${col1}"
        test_makefiledir[unitTests]="${col2}"
        test_monccomponents[unitTests]="${col3}"
        test_unittests[unitTests]="${col4}"
        test_testurl[unitTests]="${col5}"
        ((unitTests++))
    done < $file
}

runspecifictest() {
    MS="\"MONC_SRC="$1"\""
    TS="\"TEST_SRC="$2"\""
    FN=test-$3
    CPL=""

    testEnv="CRAYOS_VERSION"
    # If this is local AND not running on ARCHER then use mpif90 - otherwise use whats in the makefile
    if [ "$testmode" = "local" ] && [ -z "${!testEnv}" ] ; then CPL="FTN=mpif90" ; fi

    pushd $5 > /dev/null

    eval make $4 $MS $TS TARGET=$FN     $CPL &> tempcompile
    if [ -f $FN ];
    then
        ./$FN &> temptest
        if grep -q "SUCCESSFUL!" temptest
        then
            echo "[Result]$4:$3:$6:1" >> tempresults
            echo ">>Build "$3" tests" >> tempresults
            cat tempcompile >> tempresults
            echo ">>Execute "$3" tests" >> tempresults
            cat temptest >> tempresults
            echo "--End Message--" >> tempresults
            if [ "$testmode" = "local" ] ; then echo "[Result]$4:$3:1" ; fi
        else
            compilerPass=0
            entirePass=0
            echo "[Result]$4:$3:$6:0" >> tempresults
            echo ">>Build "$3" tests" >> tempresults
            cat tempcompile >> tempresults
            echo ">>Execute "$3" tests" >> tempresults
            cat temptest >> tempresults
            echo "--End Message--" >> tempresults
            if [ "$testmode" = "local" ] ; then echo "[Result]$4:$3:0" ; fi
        fi
    else
        compilerPass=0
        entirePass=0
        echo "[Result]$4:$3:$6:0" >> tempresults
        echo ">>Build "$3" tests" >> tempresults
        cat tempcompile >> tempresults
        echo "--End Message--" >> tempresults
        if [ "$testmode" = "local" ] ; then echo "[Result]$4:$3:0" ; fi
    fi
    rm tempcompile temptest $FN
    popd > /dev/null
    cat $5/tempresults >> $resultsFile
    rm $5/tempresults
}

main $1
