
setup() {
    # get the containing directory of this file
    # use $BATS_TEST_FILENAME instead of ${BASH_SOURCE[0]} or $0,
    # as those will point to the bats executable's location or the preprocessed file respectively
    TEST_DIR="$( cd "$( dirname "$BATS_TEST_FILENAME" )" >/dev/null 2>&1 && pwd )"
    TEST_FILES_DIR="$( cd "$( dirname "$TEST_DIR" )" >/dev/null 2>&1 && pwd )"
    TEST_FILES_DIR="$TEST_FILES_DIR/Examples"
    TEST_OUTPUT="$TEST_DIR/test_output/$BATS_TEST_NAME"
    mkdir -p $TEST_OUTPUT
    # make executables in src/ visible to PATH
    export TEST_OUTPUT
    export TEST_FILES_DIR
    PATH="$TEST_DIR/../build:$PATH"
    
    [ ${BATS_TEST_NUMBER} = 1 ] &&  echo "# Test " > ./bats.log
   
    

    INDEX=$((${BATS_TEST_NUMBER} - 1))
    echo "#####" >> ./bats.log
    echo "BATS_TEST_NAME:        ${BATS_TEST_NAME}" >> ./bats.log
    # echo "BATS_TEST_FILENAME:    ${BATS_TEST_FILENAME}" >> ./bats.log
    # echo "BATS_TEST_DIRNAME:     ${BATS_TEST_DIRNAME}" >> ./bats.log
    echo "BATS_TEST_NAMES:       ${BATS_TEST_NAMES[$INDEX]}" >> ./bats.log
    echo "BATS_TEST_DESCRIPTION: ${BATS_TEST_DESCRIPTION}" >> ./bats.log
    echo "BATS_TEST_NUMBER:      ${BATS_TEST_NUMBER}" >> ./bats.log
    # echo "BATS_TMPDIR:           ${BATS_TMPDIR}" >> ./bats.log

    echo "====" >> ./bats.log
}
teardown() {
    echo "TEST_CMD:      ${BATS_RUN_COMMAND}" >> ./bats.log

    echo -e "##### ${BATS_TEST_NAME} #####\n" >> ./bats.log
}