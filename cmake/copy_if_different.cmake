function(COPY_IF_DIFFERENT to_dir files targets tags)
   # macro to implement copy_if_different for a list of files
   # arguments - 
   #   to_dir   - this is the destination directory
   #   files    - names of the files to copy 
   #              todo: add globing. 
   #   targets  - list of targets
   #   tags     - since only the file name is used
   #              to generate rules, pre-pend a user 
   #              supplied tag to prevent duplicate rule errors. 
   set(addtargets "")
   foreach(src ${files})
       #srcfile is just the name
       get_filename_component(srcfile ${src} name) 

       set(to   ${to_dir}/${srcfile})

       add_custom_command(
           output  ${to}
           command ${cmake_command}
           args    -e copy_if_different ${src} ${to}
           comment "copying ${srcfile} ${to_dir}"
           )
       set(addtargets ${addtargets} ${to})
   endforeach()
   set(${targets} ${addtargets})
endfunction()
