message("Analysing ${BINARY_FILE} ...")

execute_process(COMMAND ${SOATL_OBJDUMP} -D ${BINARY_FILE} OUTPUT_FILE ${BINARY_FILE}.asm OUTPUT_QUIET ERROR_QUIET)

file(STRINGS ${BINARY_FILE}.asm BINARRY_ASSEMBLY)

set(INSLIST vmov vmovapd vmovaps vmovupd vmovups sqrt vsqrtpd vrsqrtpd vsqrtps vrsqrtps vsqrtsd vrsqrtsd vsqrtss vrsqrtss)

set(mova vmovapd vmovaps)
set(movu vmovupd vmovups)
set(sqrtp vsqrtpd vrsqrtpd vsqrtps vrsqrtps)
set(sqrts vsqrtsd vrsqrtsd vsqrtss vrsqrtss)
set(INSSUM mova movu sqrtp sqrts)

set(INSREPORT vmov mova movu sqrt sqrtp sqrts)

foreach(ins ${INSLIST})
  set(${ins} 0)
endforeach()

foreach(line ${BINARRY_ASSEMBLY})
  foreach(ins ${INSLIST})
    if(${line} MATCHES ".*${ins}.*")
      math(EXPR ${ins} ${${ins}}+1)
    endif()
  endforeach()
endforeach()

foreach(tosum ${INSSUM})
  set(tmplist ${${tosum}})
  set(${tosum} 0)
  foreach(ins ${tmplist})
    math(EXPR ${tosum} ${${tosum}}+${${ins}})
  endforeach()
endforeach()

foreach(ins ${INSREPORT})
  message("${ins} ${${ins}}")
endforeach()


