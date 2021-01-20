;+
; Name:
;   which
;
; Purpose:
;   Find the location of the named IDL program file.  This is
;   similar to the unix 'which' command.  
;
; Procedure:
;   The !path variable is searched for any occurrences of a file 
;   named {name}.pro 
;
;   If no files are found, search for intrinsic procedures/functions
;   with the input name and explain if found.
;
; Calling Sequence:
;   which, [ name, files=, /dlm, /show, /help, /silent ]
;
; OPTIONAL INPUT:
;   name - Character string giving the name of the IDL procedure or 
;       function.  Can also be a glob pattern, such as 'name*'.  
;       If omitted, the program will prompt
;	/dlm: Search !DLM_PATH for dynamically loaded modules instead
;       of !PATH for .pro files
;   /show: if set, will display the first file found
;   /help: print syntax
;   /silent: Don't print the locations of the files.
;
; OPTIONAL OUTPUT KEYWORDS:
;   files: An array containing the filenames
;
; EXAMPLE:
;   Find out where the procedure CURVEFIT lives.
;
;       IDL> which, 'curvefit'
;       /home/products/Linux/idl/v7_0/lib/curvefit.pro
;
;   In this case there are multiple versions on the path:
;
;       IDL> which,'mrdfits
;       Using:
;           /home/products/Linux/idlutils/v5_4_21/goddard/pro/fits/mrdfits.pro
;       Also found:
;           /home/products/Linux/idl/v7_0/lib/pro/mrdfits.pro
;   
;   Globs also work:
;       IDL> which,'*fits*'
;
; REVISION HISTORY:
;   29-MAY-94  Modified from getpro.pro by E. Deutsch
;   14-JUL-95  Fixed for IDL 4.0
;   03-DEC-2000 Added files and show keywords. Erin Sheldon
;   21-JUN-2004 Use FILE_WHICH procedure for IDL > 5.3  for significant
;       speed increase. Fixed intrinsic procedure searching. E.S.
;       Look also for dynamically loaded modules if /dlm is set
;   27-May-2009: Due to bug in idl 7.0 file_which, moved over to using 
;       straight file_search.  This now makes the dividing line 5.5 between
;       the two methods.  Erin Sheldon, BNL
;   2010-12-15: Simplify: Just check every directory in !path for the filename.
;       Dropped VMS support.  ESS, BNL
;
;   2012-4-1: made it search the current working directory first.
;        Matthew Becker, UChicago
;-

pro _which_display_results, name, files, funcfound, profound, isglob=isglob, show=show

    pname = "'"+name+"'"

    n_intrinsic = n_elements(funcfound) + n_elements(profound)
    if files[0] eq '' and n_intrinsic eq 0 then begin
        if isglob then begin
            print,'No files or intrinsic procedures matched pattern '+pname
        endif else begin
            print,'No files or intrinsic procedures matched '+pname
        endelse
    endif else begin
        _which_display_files, name, files, isglob=isglob, show=show
        _which_display_intrinsic, name, funcfound, profound
    endelse

end

pro _which_display_files, name, files, isglob=isglob, show=show
    if files[0] eq '' then return
    
    nf = n_elements(files)
    if nf eq 1 then begin
        print,files
    endif else begin

        pname = "'"+name+"'"

        for i=0L, nf-1 do begin
            if i eq 0 then begin
                if isglob then print,'Files matching pattern '+pname+':' else print,'Using:'
                print,'    '+files[i]
            endif else begin
                if i eq 1 and not isglob then print,'Also found:'
                print,'    '+files[i]
            endelse
        endfor
    endelse

    if keyword_set(show) then spawn,'more '+files[0]
end

pro _which_display_intrinsic, name, funcfound, profound
    nf = n_elements(funcfound)
    np = n_elements(profound)

    pname = "'"+name+"'"

    if nf gt 0 then begin
        print,'Intrinsic IDL functions matching '+pname+':'
        for i=0L, nf-1 do begin
            print,'    '+funcfound[i]
        endfor
    endif
    if np gt 0 then begin
        print,'Intrinsic IDL procedures matching '+pname+':'
        for i=0L, np-1 do begin
            print,'    '+profound[i]
        endfor
    endif

end


function _which_extract_globbed, globbed_files, count=count
    for i=0L, n_elements(globbed_files)-1 do begin
        f = file_search(globbed_files[i])
        if f[0] ne '' then begin
            add_arrval, f, files
        endif
    endfor
    count = n_elements(files)
    if count eq 0 then files=''
    return, files
end

function _which_find_file, name, dlm=dlm, isglob=isglob

    ; Get current IDL path of directories
    if keyword_set(dlm) then begin 
        path_string = !dlm_path
        ext = '.dlm'
    endif else begin 
        cd,c=c
        path_string = c+':'+!path
        ext = '.pro'
    endelse                 
    sep = ':'

    files = strsplit(path_string, sep, /extract)
    files = files + path_sep() + name+ext

    w = where( file_test(files), nfound ) 
    if nfound gt 0 then begin
        files=files[w]
        if isglob then begin
            files = _which_extract_globbed(files, count=nfound)
        endif
    endif else begin
        files=''
    endelse

    return, files
end

pro _which_find_intrinsic, name, funcfound, profound, isglob=isglob

    common which_intrinsic_block, funcnames, pronames
    if n_elements(funcnames) eq 0 then begin 
        proNames  = routine_info(/system)
        funcNames = routine_info(/system,/func)
    endif 

    uname = strupcase(name)

    if isglob then begin
        wfunc = where(strmatch(funcnames,uname), fcount)
        wpro = where(strmatch(pronames,uname), pcount)
    endif else begin
        wfunc = where(funcnames eq uname, fcount)
        wpro  = where(pronames  eq uname, pcount)
    endelse

    if fcount gt 0 then begin
        funcfound = funcnames[wfunc]
    endif
    if pcount gt 0 then begin
        profound = pronames[wpro]
    endif

end

pro which,proc_name,files=files,show=show,help=help, dlm=dlm, silent=silent

    if keyword_set(help) then begin
        print,'-syntax: which, proc_name, files=files, /dlm, /show, /help, /silent'
        print,'name can be a glob, e.g. "str*"'
        return
    endif 


    ; prompt for procedure name?
    if (n_params() eq 0) then begin 	     
        proc_name = ' ' 
        read,'Enter name of procedure to look for: ',proc_name     
    endif else begin
        zparcheck, 'which', proc_name, 1, 7, 0, 'Procedure name'
    endelse

    ; Don't want file extensions
    fdecomp, proc_name, disk, dir, name      
    name = strtrim( name, 2 )  

    isglob=0
    if strpos(name,'*') ne -1 then isglob=1

    files = _which_find_file(name, dlm=dlm, isglob=isglob)
    if files[0] eq '' or isglob then begin
        _which_find_intrinsic, name, funcfound, profound, isglob=isglob
    endif

    if not keyword_set(silent) then begin
        _which_display_results, name, files, funcfound, profound, isglob=isglob, show=show
    endif

end 
  
