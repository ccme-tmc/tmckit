macro(add_test_simple testname shell)
  enable_testing()
  if (shell)
    set(testname "${shell} ${testname}")
  endif (shell)
  add_test(${testname}
    ${CMAKE_COMMAND}
    -Ddir_input=${CMAKE_CURRENT_SOURCE_DIR}
    -Dname=${testname}
    -Dcmd=${CMAKE_CURRENT_SOURCE_DIR}/${testname}
    -Dfile_input=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.in
    -P ${CMAKE_SOURCE_DIR}/cmake/runtest.cmake
    )
endmacro(add_test_simple)


macro(add_test_text testname shell)
  enable_testing()
  if (shell)
    set(testname "${shell} ${testname}")
  endif (shell)
  add_test(${testname}
    ${CMAKE_COMMAND}
    -Ddir_input=${CMAKE_CURRENT_SOURCE_DIR}
    -Dname=${testname}
    -Dcmd=${CMAKE_CURRENT_SOURCE_DIR}/${testname}
    -Dfile_input=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.in
    -Dfile_ref=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.out
    -P ${CMAKE_SOURCE_DIR}/cmake/runtest.cmake
    )
endmacro(add_test_text)

macro(add_test_floatdata testname shell)
  enable_testing()
  if (shell)
    set(testname "${shell} ${testname}")
  endif (shell)
  add_test(${testname}
    ${CMAKE_COMMAND}
    -Ddir_input=${CMAKE_CURRENT_SOURCE_DIR}
    -Dname=${testname}
    -Dcmd=${CMAKE_CURRENT_SOURCE_DIR}/${testname}
    -Dfile_input=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.in
    -Dfile_ref=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.out
    -Dcmd_comp=${CMAKE_SOURCE_DIR}/cmake/comp_num.py
    -P ${CMAKE_SOURCE_DIR}/cmake/runtest.cmake
    )
endmacro(add_test_floatdata)

macro(add_test_zero testname shell)
  enable_testing()
  if (shell)
    set(testname "${shell} ${testname}")
  endif (shell)
  add_test(${testname}
    ${CMAKE_COMMAND}
    -Ddir_input=${CMAKE_CURRENT_SOURCE_DIR}
    -Dname=${testname}
    -Dcmd=${CMAKE_CURRENT_SOURCE_DIR}/${testname}
    -Dfile_input=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.in
    -Dcmd_comp=${CMAKE_SOURCE_DIR}/cmake/comp_num.py
    -P ${CMAKE_SOURCE_DIR}/cmake/runtest.cmake
    )
endmacro(add_test_zero)

macro(add_test_dir testname shell)
  enable_testing()
  if (shell)
    set(testname "${shell} ${testname}")
  endif (shell)
  add_test(${testname}
    ${CMAKE_COMMAND}
    -Ddir_input=${CMAKE_CURRENT_SOURCE_DIR}
    -Dname=${testname}
    -Dcmd=${CMAKE_CURRENT_SOURCE_DIR}/${testname}
    -Dfile_input=${CMAKE_CURRENT_SOURCE_DIR}/${testname}.in
    -Ddir_out=${CMAKE_CURRENT_BINARY_DIR}
    -Ddir_ref=${CMAKE_CURRENT_SOURCE_DIR}/${testname}-ref
    -Dcmd_comp=${CMAKE_SOURCE_DIR}/cmake/comp_num.py
    -P ${CMAKE_SOURCE_DIR}/cmake/runtest.cmake
    )
endmacro(add_test_dir)
