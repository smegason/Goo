#!/usr/bin/env sh

SIGINT_PIDFILE=/tmp/pytest-blender-integration-sigint.pid

_testSigIntPropagation() {
    rm -f $SIGINT_PIDFILE

    _blender_executable_arg=""
    if [ -n "$BLENDER_EXECUTABLE" ]; then
      _blender_executable_arg="--blender-executable $BLENDER_EXECUTABLE"
    fi;
    SIGINT_PIDFILE="$SIGINT_PIDFILE" \
      python3 -m pytest -svv --noconftest \
      $_blender_executable_arg \
      tests/integration/sigint.py &
    # wait some time to start the test suite execution
    sleep 3

    # send sighup
    pid="$(cat /tmp/pytest-blender-integration-sigint.pid)"
    assertNotEquals "PID not found" "" "$pid"
    kill -s $1 $pid
    sleep 1

    # SIGINT needs two signals
    if [ "$1" -eq 2 ]; then
      kill -s 2 $pid
      sleep 1
    fi;

    # check process was killed
    ps -p $pid | grep foo

    # remove pidfile
    rm -f $SIGINT_PIDFILE
}

testSigIntPropagation() {
  _testSigIntPropagation 2
}

testSigHupPropagation() {
  _testSigIntPropagation 1
}

testSigTermPropagation() {
  _testSigIntPropagation 15
}

oneTimeTearDown() {
  rm -f $SIGINT_PIDFILE
}

prepare() {
  if [ ! -f "shunit2" ]; then
    curl -sSL https://raw.githubusercontent.com/kward/shunit2/master/shunit2 \
      -o shunit2
  fi
}

prepare && . ./shunit2