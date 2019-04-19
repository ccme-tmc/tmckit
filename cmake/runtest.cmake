set (file_testout ${CMAKE_CURRENT_BINARY_DIR}/${name}.out)
set (file_testerr ${CMAKE_CURRENT_BINARY_DIR}/${name}.err)

#Choose compare method
#Note: if one would like to compare folders, then stdout/stderr files will not be used
if (cmd_comp)
    if (dir_ref)
        set (cmd_compare "${cmd_comp}" -r "${dir_out}" "${dir_ref}")
    else()
        set (cmd_compare "${cmd_comp}" "${file_testout}" "${file_ref}")
    endif(dir_ref)
else()
    set (cmd_compare "${CMAKE_COMMAND}" -E compare_files "${file_testout}" "${file_ref}")
endif(cmd_comp)

if (WIN32)
  set(cmd cmd.exe /c ${cmd})
  set(cmd_compare cmd.exe /c ${cmd_compare})
endif()

if (file_input AND EXISTS "${file_input}")
  message("Test command: ${cmd} ${file_input} > ${file_testout}")
  execute_process(
    COMMAND ${cmd} ${dir_input}
    RESULT_VARIABLE result 
    INPUT_FILE ${file_input}
    OUTPUT_FILE ${file_testout}
    ERROR_FILE ${file_testerr}
    TIMEOUT 1000
    )
else()
  message("Test command: ${cmd} > ${file_testout}")
  execute_process(
    COMMAND ${cmd} ${dir_input}
    RESULT_VARIABLE result 
    OUTPUT_FILE ${file_testout}
    ERROR_FILE ${file_testerr}
    TIMEOUT 1000
    )
endif()

if(result) #Only empty means success
  message(SEND_ERROR "error in test '${name}'; command: '${cmd}'; output: ${result}")
endif()

if(file_ref STREQUAL "" AND dir_ref STREQUAL "")
    message("Compare commad1: ${cmd_compare}")

    execute_process(
      COMMAND ${cmd_compare}
# COMMAND cmd.exe /c "${CMAKE_COMMAND}" -E compare_files "${file_testout}" "${file_ref}" 
      RESULT_VARIABLE result 
      OUTPUT_VARIABLE out
      ERROR_VARIABLE err
      TIMEOUT 1000
      )

    if(result) #Only empty means success
      message(SEND_ERROR "output does not match in test '${name}': ${err}; command ${cmd}  : shell return: ${result} : output : ${out}")
    endif()
endif()

