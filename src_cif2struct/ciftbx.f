!
!
!    \ | /            /##|    @@@@  @   @@@@@   |      |              @@@
!     \|/ STAR       /###|   @      @   @     __|__    |             @   @
!  ----*----        /####|  @       @   @@@@    |      |___  __  __     @
!     /|\          /#####|   @      @   @       |      |   \   \/      @
!    / | \         |#####|    @@@@  @   @       \___/  \___/ __/\__  @@@@@
!                  |#####|________________________________________________
!                 ||#####|                 ___________________            |
!        __/|_____||#####|________________|&&&&&&&&&&&&&&&&&&&||          |
!<\\\\\\\\_ |_____________________________|&&& 16 Jun 1998 &&&||          |
!          \|     ||#####|________________|&&&&&&&&&&&&&&&&&&&||__________|
!                  |#####|
!                  |#####|                Version 2.6.2 Release
!                  |#####|
!                 /#######\ 
!                |#########|
!                    ====
!                     ||
!           An extended tool box of fortran routines for manipulating CIF data.
!                     ||
!                     ||  CIFtbx Version 2
!                     ||        by
!                     ||
!                     ||  Sydney R. Hall (syd@crystal.uwa.edu.au)
!                     ||  Crystallography Centre
!                     ||  University of Western Australia
!                     ||  Nedlands 6009, AUSTRALIA
!                     ||
!                     ||       and
!                     ||
!                     ||  Herbert J. Bernstein (yaya@bernstein-plus-sons.com)
!                     ||  Bernstein + Sons
!                     ||  5 Brewster Lane
!                     ||  Bellport, NY 11713, U.S.A.
!                     ||
! The latest program source and information is available from:
!                     ||
! Em: syd@crystal.uwa.edu.au       ,-_|\      Sydney R. Hall
! sendcif@crystal.uwa.edu.au      /     \     Crystallography Centre
! Fx: +61 9 380 1118  ||      --> *_,-._/     University of Western Australia
! Ph: +61 9 380 2725  ||               v      Nedlands 6009, AUSTRALIA
!                     ||
!                     ||
!_____________________||_____________________________________________________
!
! This is a version of CIFtbx which has been extended to work with DDL 2
! and mmCIF as well as with DDL 1.4 and core CIF dictionaries.  CIFtbx
! version 1 was written by Sydney R. Hall (see Hall, S. R., "CIF Applications
! IV.  CIFtbx: a Tool Box for Manipulating CIFs,"  J. Appl. Cryst (1993). 26,
! 482-494.  The revisions for version 2 were done by Herbert J. Bernstein
! and Sydney R. Hall (see Hall, S. R. and Bernstein, H. J., "CIFtbx 2:
! Extended Tool Box for Manipulating CIFs," J. Appl. Cryst.(1996). 29,
! 598-603) 
!
!___________________________________________________________________________
!
!
!    GENERAL TOOLS
!
!
!    init_      Sets the device numbers of files.   (optional)
!               [logical function always returned .true.]
!
!               <input CIF dev number> Set input CIF device     (def=1)
!
!               <output CIF dev number>Set output CIF device    (def=2)
!
!               <diracc dev number>    Set direct access formatted
!                                      scratch device number    (def=3)
!
!               <error  dev number>    Set error message device (def=6)
!
!
!
!    dict_      Requests a CIF dictionary be used for various data checks.
!               [logical function returned as .true. if the name dictionary
!               was opened; and if the check codes are recognisable.  The
!               data item names used in the first dictionary loaded are
!               considered to be preferred by the user to aliases found
!               in dictionaries loaded in later calls.  On exit from dict_
!               the variable dicname_ is either equal to the filename, or,
!               if the dictionary had a value for the tag dictionary_name
!               of dictionary.title, dicname_ is set to that value.
!               The variable dicver_ is blank or set to the value of
!               _dictionary_version or of _dictionary.version  The check codes
!               'catck' and 'catno' turn on and off checking of dictionary
!               catgeory conventions.  The default is 'catck'.  Three check
!               codes control the handling of tags from the current dictionary
!               which duplicate tags from a dictionary loaded earlier.  These
!               codes ('first', 'final' and 'nodup') have effect only for the
!               current call to dict_  The default is 'first'.]
!
!               <dictionary filename>  A CIF dictionary in DDL format
!                                      or blank if just setting flags
!                                      or resetting the dictionary
!
!               <check code string>    The codes specifying the types of 
!                                      checks to be applied to the CIF.
!
!                                      'valid'  data name validation check.
!                                      'dtype'  data item data type check.
!                                      'catck'  check datanames against
!                                               categories
!                                      'catno'  don't check datanames against
!                                               categories
!                                      'first'  accept first dictionary's 
!                                               definitions of duplicate tags
!                                      'final'  accept final dictionary's
!                                               definitions of duplicate tags
!                                      'nodup'  do not accept duplicate tag
!                                               definitions
!                                      'reset'  switch off checking flags
!                                      'close'  close existing dictionaries
!
!___________________________________________________________________________
!
!
!   CIF ACCESS TOOLS  ("the get_ing commands")
!
!
!
!    ocif_      Opens the CIF containing the required data.
!               [logical function returned .true. if CIF opened]
!
!               <CIF filename>        A blank name signals that the
!                                     currently open input CIF file
!                                     will be read.
!
!
!
!    data_      Identifies the data block containing the data to be requested. 
!               [logical function returned .true. if block found]
!
!               <data block name>     A blank name signals that the next
!                                     encountered block is used (the block
!                                     name is stored in the variable bloc_).
!
!
!    bkmrk_     Saves or restores the current position so that data from 
!               elsewhere in the cif can be examined.
!               [logical function returned as .true. on save if there was
!               room in internal storage to hold the current position, .true.
!               on restore if the bookmark number used was valid.  If the
!               argument is zero, the call is to save the position and return
!               the bookmark number in the argument.  If the argument is
!               non-zero, the call is to restore the position saved for the
!               bookmark number given.  The bookmark and the argument are
!               cleared.  The position set on return allow reprocessing of
!               the data item or loop row last processed when the bookmark
!               was placed.
!
!               NOTE:  All bookmarks are cleared by a call to data_]
!
!               <integer variable>    Bookmark number
!
!
!    find_      Find the location of the requested item in the CIF.
!               [The argument "name" may be a data item name, blank
!               for the next such item.  The argument "type" may be
!               blank for unrestricted acceptance of any non-comment
!               string (use cmnt_ to see comments), including loop headers,
!               "name" to accept only the name itself and "valu"
!               to accept only the value, or "head" to position to the
!               head of the CIF.  Except when the "head" is requested,
!               the position is left after the data item provided.  If the
!               item found is of type "name", posnam_ is set, otherwise, 
!               posval_]
!
!               <data item name>      A blank name signals that the next
!                                     item of the type specified is needed
!
!               <data item type>      blank, 'head', 'name' or 'valu'
!
!               <character variable>  Returned string is of length long_.
!
!
!
!    test_      Identify the data attributes of the named data item.
!               [logical function returned as .true. if the item is present or
!               .false. if it is not. The data attributes are stored in the
!               common variables list_, type_, dictype_, diccat_ and dicname_. 
!               The values in dictype_, diccat_ and dicname_ are valid
!               whether or not the data item is found in the input CIF, as
!               long as the named data item is found in the dictionaries
!               declared by calls to dict_.  The data item name found
!               in the input CIF is stored in tagname_.  The appropriate
!               column numbers are stored in posnam_, posval_, posend_ and (for
!               numbers) in posdec_.  The quotation mark, if any, used is
!               stored in quote_.
!
!               list_ is an integer variable containing the sequential number
!               of the loop block in the data block. If the item is not within
!               a loop structure this value will be zero.
!
!               type_ is a character*4 variable with the possible values:
!                      'numb'  for number data
!                      'char'  for character data
!                      'text'  for text data
!                      'null'  if data missing or '?' or '.'
!                              also used for blank quoted fields if
!                              nblank_ is true
!
!               dictype_ is a character*(NUMCHAR) variable with the type code
!               given in the dictionary entry for the named data item.  If
!               no dictionary was used, or no type code was specified, this
!               field will simply agree with type_.  If a dictionary was used,
!               this type may be more specific than the one given by type_.
!
!               diccat_ is a character*(NUMCHAR) variable with the category
!               of the named data item, or '(none)'
!
!               dicname_ is a character*(NUMCHAR) variable with the name of
!               the data item which is found in the dictionary for the
!               named data item.  If alias_ is .true., this name may
!               differ from the name given in the call to test_.  If alias_
!               is .false. or no preferred alias is found, dicname_ agrees with
!               the data item name. 
!
!               tagname_ is a character*(NUMCHAR) variable with the name
!               of the data item as found in the input CIF.  It will be
!               blank if the data item name requested is not found in the
!               input CIF and may differ from the data item name provided
!               by the user if the name used in the input CIF is an
!               alias of the data item name and alias_ is .true.
!
!               posnam_, posval_, posend_  and posdec_ are integer variables
!               which may be examined if information about the horizontal
!               position of the name and data read are needed.  posnam_ is the
!               starting column of the data name found (most often 1).
!               posval_ is the starting column of the data value.  If the
!               field is numeric, then posdec_ will contain the effective
!               column number of the decimal point.  For whole numbers, the
!               effective position of the decimal point is one column to the
!               right of the field.  posend_ contains the ending column of the
!               data value.
!
!               quote_ is a character*1 varibale which may be examined to
!               determine if a quotation character was used on character data.]
!
!               <data name>           Name of the data item to be tested.
!
!
!
!    name_      Get the NEXT data name in the current data block.
!               [logical function returned as .true. if a new data name exists
!               in the current data block, and .false. when the end of the data
!               block is reached.]
!
!               <data name>           Returned name of next data item in block.
!
!
!
!    numb_      Extracts the number and its standard deviation (if appended).
!               [logical function returned as .true. if number present. If
!               .false. arguments 2 and 3 are unaltered. If the esd is not
!               attached to the number argument 3 is unaltered.]
!
!               <data name>           Name of the number sought.
!
!               <real variable>       Returned number.
!
!               <real variable>       Returned standard deviation.
!
!
!
!    numd_      Extracts the number and its standard deviation (if appended)
!               as double precision variables.
!               [logical function returned as .true. if number present. If
!               .false. arguments 2 and 3 are unaltered. If the esd is not
!               attached to the number argument 3 is unaltered.]
!
!               <data name>           Name of the number sought.
!
!               <double precision variable>
!                                     Returned number.
!
!               <double precision variable>
!                                     Returned standard deviation.
!
!
!
!    char_      Extracts character and text strings.
!               [logical function returned as .true. if the string is present.
!               Note that if the character string is text this function is 
!               called repeatedly until the logical variable text_ is .false.
!               Non-text blank (quoted blanks) or empty ('' or "") fields
!               are converted by char to a null field, if nblank_ is true.]
!
!               <data name>           Name of the string sought.
!
!               <character variable>  Returned string is of length long_.
!
!
!    cmnt_      Extracts the next comment from the data block.
!               [logical function returned as .true. if a comment is present.
!               The initial comment character "#" is _not_ included in the
!               returned string.  A completely blank line is treated as
!               a comment.]
!
!               <character variable>  Returned string is of length long_.
!
!
!
!    purge_     Closes existing data files and clears tables and pointers.
!               [subroutine call]        
!
!____________________________________________________________________________
!
!
!
!   CIF CREATION TOOLS ("the put_ing commands")
!
!
!
!    pfile_     Create a file with the specified file name.
!               [logical function returned as .true. if the file is opened.
!               The value will be .false. if the file already exists.]
!
!               <file name>           Blank for use of currently open file
!
!
!
!    pdata_     Put a data block command into the created CIF. 
!               [logical function returned as .true. if the block is created.
!               The value will be .false. if the block name already exists.
!               Produces a save frame instead of a data block if the
!               variable saveo_ is true during the call.  No block duplicate
!               check is made for a save frame.]
!
!               <block name>
!
!
!
!    ploop_     Put a loop_ data name into the created CIF.             
!               [logical function returned as .true. if the invocation 
!               conforms with the CIF logical structure.  If pposval_ 
!               is non-zero, the "loop_" header is positioned to 
!               that column.  If pposnam_ is non-zero, the data name is 
!               positioned to that column.]
!
!               <data name>         If the name is blank on the first call
!                                   of a loop, only the "loop_" is placed.
!
!
!
!    pchar_     Put a character string into the created CIF.             
!               [logical function returned as .true. if the name is unique,
!               AND, if dict_ is invoked, is a name defined in the dictionary, 
!               AND, if the invocation conforms to the CIF logical structure.
!               The action of pchar_ is modified by the variables pquote_ and
!               nblanko_.  If pquote_ is non-blank, it is used as a quotation
!               character for the string written by pchar_.  The valid values
!               are '''', '"', and ';'.  In the last case a text field is
!               written.  If the string contains a matching character to the
!               value of quote_, or if quote_ is not one of the valid
!               quotation characters, a valid, non-conflicting quotation
!               character is used.  Except when writing a text field, if 
!               nblanko_ is true, pchar_ converts a blank string to 
!               an unquoted period.]
!
!               <data name>         If the name is blank, do not output name.
!
!               <character string>  A character string of MAXBUF chars or less.
!
!
!
!    pcmnt_     Puts a comment into the created CIF.
!               [logical function returned as .true.  The comment character
!               "#" should not be included in the string.  A blank comment
!               is presented as a blank line without the leading "#"].
!
!               <character string>  A character string of MAXBUF chars or less.
!
!
!    pnumb_     Put a single precision number and its esd into the created CIF.
!               [logical function returned as .true. if the name is unique,
!               AND, if dict_ is invoked, is a name defined in the dictionary, 
!               AND, if the invocation conforms to the CIF logical structure.
!               The number of esd digits is controlled by the variable
!               esdlim_]
!
!               <data name>         If the name is blank, do not output name.
!
!               <real variable>     Number to be inserted.
!
!               <real variable>     Esd number to be appended in parentheses.
!
!
!    pnumd_     Put a double precision number and its esd into the created CIF.
!               [logical function returned as .true. if the name is unique,
!               AND, if dict_ is invoked, is a name defined in the dictionary, 
!               AND, if the invocation conforms to the CIF logical structure.
!               The number of esd digits is controlled by the variable
!               esdlim_]
!
!               <data name>         If the name is blank, do not output name.
!
!               <double precision variable>  
!                                   Number to be inserted.
!
!               <double precision variable>  
!                                   Esd number to be appended in parentheses.
!
!
!
!    ptext_     Put a character string into the created CIF.             
!               [logical function returned as .true. if the name is unique,
!               AND, if dict_ is invoked, is a name defined in the dictionary, 
!               AND, if the invocation conforms to the CIF logical structure.]
!               ptext_ is invoked repeatedly until the text is finished. Only
!               the first invocation will insert a data name.
!
!               <data name>         If the name is blank, do not output name.
!
!               <character string>  A character string of MAXBUF chars or less.
!
!
!    prefx_     Puts a prefix onto subsequent lines of the created CIF.
!               [logical function returned as .true.  The second argument
!               may be zero to suppress a previously used prefix, or
!               greater than the non-blank length of the string to force
!               a left margin.  Any change in the length of the prefix 
!               string flushes pending partial output lines, but does _not_
!               force completion of pending text blocks or loops.
!               This function allows the CIF output functions to be used 
!               within what appear to be text fields to support annotation 
!               of a CIF. ]
!
!               <character string>  A character string of MAXBUF chars or less.
!
!               <integer variable>  The length of the prefix string to use.
!
!
!
!
!    close_     Close the creation CIF. MUST be used if pfile_ is used.
!               [subroutine call]
!
!
!____________________________________________________________________________
!
!
!
!....The CIF tool box also provides variables for data access control:
! 
!
!    alias_      Logical variable: if left .true. then all calls to
!                CIFtbx functions may use aliases of data item names.
!                The preferred synonym from the dictionary will be
!                subsituted internally, provided aliased data names were
!                supplied by an input dictionary (via dict_).  The
!                default is .true., but alias_ may be set to .false.
!                in an application.
!
!    aliaso_     Logical variable: if set .true. then cif output 
!                routines will convert aliases to the names to preferred
!                synonyms from the dictionary.  The default is .false., but
!                aliaso_ may be set to .true. in an application.  The
!                setting of aliaso_ is independent of the setting of
!                alias_.
!
!    align_      Logical variable signals alignment of loop_ lists during
!                the creation of a CIF. The default is .true.
!
!    append_     Logical variable:  if set .true. each call to ocif_ will
!                append the information found to the current cif.  The default
!                is .false.
!
!    bloc_       Character*(NUMCHAR) variable: the current block name.
!
!    decp_       Logical variable: set when processing numeric input, .true.
!                if there is a decimal point in the numeric value, .false.
!                otherwise
!
!    dictype_    Character*(NUMCHAR) variable: the precise data type code
!                (see test_)
!
!    diccat_     Character*(NUMCHAR) variable: the category (see test_)
!
!    dicname_    Character*(NUMCHAR) variable: the root alias (see test_) of
!                name of the dictionary just loaded (see dict_)
!
!    dicver_     Character*(NUMCHAR) variable: the version of the dictionary
!                just loaded (see dict_)
!
!    esdlim_     Integer variable:  Specifies the upper limit of esd's
!                produced by pnumb_, and, implicitly, the lower limit.
!                The default value is 19, which limits esd's to the range
!                2-19.  Typical values of esdlim_ might be 9 (limiting
!                esd's to the range 1-9), 19, or 29 (limiting esd's
!                to the range 3-29).  If esdlim_ is given as a negative
!                value, the upper limit of esd's is the absolute value
!                of esdlim_ and the lower limit is 1.
!
!    esddig_     Integer variable:  The number of esd digits in the last
!                number read from a CIF.  Will be zero if no esd
!                was given.
!
!    file_       Character*(MAXBUF) variable: the filename of the current file.
!
!    glob_       Logical variable signals that the current data block
!                is actually a global block (.true. for a global block).
!
!    globo_      Logical variable signals that the output data block from
!                pdata_ is actually a global block (.true. for a global block).
!
!    line_       Integer variable: Specifies the input/output line limit
!                for processing a CIF. The default value is 80 characters.
!                This may be set by the program. The max value is MAXBUF
!                which has a default value of 200.
!
!    list_       Integer variable: the loop block number (see test_).
!
!    long_       Integer variable: the length of the data string in strg_.
!
!    longf_      Integer variable: the length of the filename in file_.
!
!    loop_       Logical variable signals if another loop packet is present.
!
!    lzero_      Logical variable: set when processing numeric input, .true.
!                if the numeric value is of the form [sign]0.nnnn rather than
!                [sign].nnnn, .false. otherwise 
!
!    nblank_     Logical variable: if set .true. then all calls to
!                to char_ or test_ which encounter a non-text quoted blank 
!                will return the type as 'null' rather than 'char'.
!
!    nblanko_    Logical variable: if set .true. then cif output 
!                routines will convert quoted blank strings to an
!                unquoted period (i.e. to a data item of type null).
!
!    pdecp_      Logical variable: if set .true. then cif numeric output
!                routines will insert a decimal point in all numbers written by
!                pnumb_ or pnumbd_.  If set .false. then a decimal point will be
!                written only when needed.  The default is .false.
!
!    pesddig_    Integer variable: if set non-zero, and esdlim_ is negative,
!                controls the number of digits for esd's produced by
!                pnumb_ and pnumd_
! 
!    plzero_     Logical variable: if set .true. then cif numeric output
!                routines will insert a zero before a leading decimal point,
!                The default is .false.
!
!    pposdec_    Integer variable giving the position of the decimal point
!                for the next number to be written.  This acts very much like
!                a decimal centered tab in a word processor, to help align
!                columns of number on a decimal point, if a decimal point
!                is present.
!
!    pposend_    Integer variable giving the ending column of the next
!                number or quoted character value to be written.  Used to
!                pad with zeros or blanks.
!
!    pposnam_    Integer variable giving the starting column of the next
!                name or comment or data block to be written.
!
!    pposval_    Integer variable giving the starting column of the next
!                data value to be written by pchar_, pnumb_ or pnumd_.
!                Also used to set the position of the initial "loop_"
!                in a ploop_ call or to set the position of a terminal "save_"
!                for a save frame in a pdata_ call for which saveo_ is .true.
!
!    posdec_     Integer variable giving the position of the decimal point
!                for the last number read, if a decimal point was present.
!
!    posend_     Integer variable giving the ending column of the last
!                data value read, not including a terminal quote.
!
!    posnam_     Integer variable giving the starting column of the last
!                name or comment or data block read.
!
!    posval_     Integer variable giving the starting column of the last
!                data value read.  Also reports the column of the
!                terminal "save_" of a save frame.
!
!    pquote_     Character variable giving the quotation symbol to be
!                used for the next string written.
!
!    precn_      Integer variable:  Reports the record number of the last
!                line written to the output cif.  Set to zero by init_.  Also
!                set to zero by pfile_ and close_ if the output cif file name
!                was not blank.
!
!    ptabx_      Logical variable signals tab character expansion to blanks 
!                during the creation of a CIF. The default is .true.
!
!    quote_      Character variable giving the quotation symbol found
!                delimiting the last string read.
!
!    recbeg_     Integer variable:  Gives the record number of the first
!                record to be used.  May be changed by the user to restrict
!                access to a CIF.
!
!    recend_     Integer variable:  Gives the record number of the last
!                record to be used.  May be changed by the user to restrict
!                access to a CIF.
!
!    recn_       Integer variable:  Reports the record number of the last
!                line read from the direct access copy of the input cif.
!
!    save_       Logical variable signals that the current data block
!                is actually a save-frame (.true. for a save-frame).
!
!    saveo_      Logical variable signals that the output data block from
!                pdata_ is actually a save-frame (.true. for a save-frame).
!
!    strg_       Character*(MAXBUF) variable: the current data item.
!
!    tabl_       Logical variable signals tab-stop alignment of output 
!                during the creation of a CIF. The default is .true.
!
!    tabx_       Logical variable signals tab character expansion to blanks 
!                during the reading of a CIF. The default is .true.
!
!    tbxver_     Character*32 variable: the CIFtbx version and date
!                in the form 'CIFtbx version N.N.N, DD MMM YY '
!
!    text_       Logical variable signals if another text line is present.
!
!    type_       Character*4 variable: the data type code (see test_).
!
!
!
!_____________________________________________________________________________
!
!
! >>>>>> Set the device numbers.
!
         function init_(devcif,devout,devdir,deverr)
!
         logical   init_
         include   'ciftbx.sys'
         integer   devcif,devout,devdir,deverr
         integer   ii,kdig
         real      ytest
         double precision ztest         
!
         init_=.true.
         cifdev=devcif
         outdev=devout
         dirdev=devdir
         errdev=deverr
         recn_=0
         precn_=0
!
!        recompute decimal single precision precision
!        This is found by computing the smallest power of
!        10 which, when added to 1, produces a change
!        and then backing off by 1
!
         decprc = .1
         do ii = 1,6
         ytest = 1.+decprc/10.
         if (ytest.eq.1.) go to 100
         decprc = decprc/10.
         enddo
100      continue
         decprc=decprc*10.
!
!        recompute decimal double precision precision
!
         kdig = 1
         dpprc = .1D0
         do ii = 1,15
         ztest = 1.D0+dpprc/10.
         if (ztest.eq.1.D0) go to 200
         dpprc = dpprc/10.D0
         kdig = kdig+1
         enddo
200      continue
         dpprc=dpprc*10.D0
         write(ndpfmt,'(5h(d30.,i2,1h))') kdig-1
!
!        recompute decimal single precision minimum power of ten
!
         decmin = .1
         do ii = 1,37
         ytest = decmin/10.
         if (ytest.eq.0.) go to 300
         decmin = decmin/10.
         enddo
300      continue
!
!        recompute decimal double precision minimum power of 10
!        and its log base 10 (minexp)
!
         dpmin = .1D0
         minexp = -1
         do ii = 1,307
         ztest = dpmin/10.
         if (ztest.eq.0.D0) go to 400
         dpmin = dpmin/10.D0
         minexp = minexp-1
         enddo
400      continue
         call clearfp
         return
         end
!
!
!
!
!
! >>>>>> Read a CIF dictionary and prepare for checks
!
         function dict_(fname,checks)
!
         logical   dict_
         logical   ocif_
         logical   data_
         logical   char_
         logical   test_
         integer   lastnb
         include  'ciftbx.sys'
         character locase*(MAXBUF)
         character fname*(*),checks*(*)
         character temp*80,codes(9)*5,name*(MAXBUF),bxname*(NUMCHAR)
         character bcname*(NUMCHAR),biname*(NUMCHAR),bname*(NUMCHAR)
         character baname*(NUMCHAR),ganame*(NUMCHAR),btname*(NUMCHAR)
         character batag*(NUMCHAR)
         character riname*(NUMCHAR),rfname*(NUMCHAR)
         character xdicnam*(NUMCHAR)
         character xdicver*(NUMCHAR)
         character*3 ovchk, otchk
         integer   nrecds,recends,recbegs
         integer   lbcname,lbaname,lbtname,lbname
         integer   lriname,lrfname
         integer   kdict,kadict,ifind,jfind,iafind,jck,ick
         integer   i,j,nmatch,mycat,ksmatch,ii,jj,idstrt,icstrt,kdup
         integer   nmycat

!
!        Control flags for matching categories, names and types
!
!        icloop is the loop number of the block for the
!        current category
!        ictype is the type of the current category
!          0 - none found yet
!          1 - _item.category.id
!          2 - _category
!          3 - _category.id
!        inloop is the loop number of the block for the
!        current name
!        intype is the type of the current name
!          0 - none found yet
!          1 - _item.name
!          2 - _name
!        ialoop is the loop number of the block for the
!        current alias
!        iatype is the type for the current alias
!          0 - none found yet
!          1 - _item_aliases.alias_name
!        itloop is the loop number of the block for the
!        current type
!        ittype is the type of the current type
!          0 - none found yet
!          1 - _item_type.code
!          2 - _type
!        iritype is the type of the current related item
!          0 - none found yet
!          1 - _item_related.related_name
!          2 - _related_item
!        irftype is the type of the current related item function
!          0 - none found yet
!          1 - _item_related.function_code
!          2 - _related_function
!

         integer icloop,ictype,inloop,intype,ialoop,iatype, &
       itloop,ittype,iriloop,iritype,irfloop,irftype,icktype
!
         character*4 map_type(12),map_to(12),mapped
         character*(NUMCHAR) dt(2),dv(2),ct(3),nt(2),at(1),tt(2)
         character*(NUMCHAR) ri(2),rf(2),ck(2)
         data map_type &
         /'floa','int ','yyyy','symo','ucha','ucod','name','idna', &
          'any ','code','line','ulin'/
         data map_to &
         /'numb','numb','char','char','char','char','char','char', &
          'char','char','char','char'/
         data ri &
            /'_item_related.related_name      ', &
             '_related_item                   '/
         data rf &
            /'_item_related.function_code     ', &
             '_related_function               '/
         data dt &
            /'_dictionary.title               ', &
             '_dictionary_name                '/
         data dv &
            /'_dictionary.version             ', &
             '_dictionary_version             '/
         data ct &
            /'_item.category_id               ', &
             '_category                       ', &
             '_category.id                    '/
         data nt &
            /'_item.name                      ', &
             '_name                           '/
         data at &
            /'_item_aliases.alias_name        '/
         data tt &
            /'_item_type.code                 ', &
             '_type                           '/
         data ck &
            /'_category_key.name              ', &
             '_list_reference                 '/

!
         data codes /'valid','dtype','reset','close', &
             'catck','catno','nodup', &
             'final','first'/
!
         nrecds=nrecd
         recbegs=recbeg_
         recends=recend_
         if(append_) then
           recbeg_=nrecd
         endif
!
!        Initialize kdup to 0 ('final')
!
         kdup = 0
!
!        initialize both xdicnam and xdicver to blank
!
         xdicnam = ' '
         xdicver = ' '
!
!        preserve entry values of tcheck and vcheck in case dict fails
!
         otchk = tcheck
         ovchk = vcheck
!
!....... Are the codes OK
!
         temp=locase(checks)
         i=0         
120      i=i+1
         if(i.ge.80)                 goto 190
         if(temp(i:i).eq.' ')        goto 120
         do 150 j=1,7
         if(temp(i:i+4).eq.codes(j)) goto 170
150      continue
         dict_=.false.
         goto 500
170      i=i+4
         if(j.eq.1) then
           vcheck='yes'
           goto 120
         endif
         if(j.eq.2) then
           tcheck='yes'
           goto 120
         endif
         if(j.eq.3) then
           vcheck = 'no '
           tcheck = 'no '
           goto 120
         endif
         if(j.eq.4) then
           vcheck = 'no '
           tcheck = 'no '
           catchk = 'yes'
           ndcname = 0
           ndict = 0
           if(nname.gt.0) then
           do 180 i = 1,nname
             dtype(i)=' '
             dxtyp(i)=' '
             cindex(i)=0
             ddict(i)=0
180        continue
           endif
           dict_=.true.
           goto 500
         endif
         if (j.eq.5) then
           catchk = 'yes'
           goto 120
         endif
         if (j.eq.6) then
           catchk = 'no '
           goto 120
         endif
         kdup=j-8
         goto 120
!
!        if no category names have been loaded, clean up
!        the hash table for dictionary category names
!
190      if(ndcname.eq.0) then
           call hash_init(dcname,dcchain,NUMDICT,ndcname,dchash, &
           NUMHASH)
         endif
         icstrt=ndcname
!
!        if no dictionary names have been loaded, clean up
!        the hash table for dictionary names
!
         if(ndict.eq.0) then
           call hash_init(dicnam,dicchain,NUMDICT,ndict,dichash, &
           NUMHASH)
         endif
         idstrt=ndict
!
!....... Open and store the dictionary
!
         dict_=.true.
         if(fname.eq.' ')            goto 500
         if(nname.gt.0) call err(' Dict_ must precede ocif_')
         dict_=ocif_(fname)
         if(.not.dict_)              goto 500
         dictfl='yes'
!
!        At this point is is proper to update xdicnam to fname
!
         xdicnam = fname         
!
!....... Loop over data blocks; extract _name's, _type etc.
!
200      if(.not.data_(' '))         goto 400
         if(bloc_(1:1).eq.'_'.or.glob_.or.bloc_.eq.' ') then
           bname=locase(bloc_)
         else
           bname='_'//locase(bloc_)
         endif
         lbname=max(1,lastnb(bname))
!
!        see if this is a dictionary defining block
!
         do i = 1,2
           if(char_(dt(i),name)) then
             xdicnam = name(1:max(1,long_))
             do j = 1,2
               if(test_(dv(j))) then
                 xdicver = strg_(1:max(1,long_))
                 goto 200
               endif  
             enddo
             goto 200
           endif
         enddo
!
!dbg     WRITE(6,*) ndict,bloc_
!
!        Analyze loop structure for categories, names and types
!
!
!        initalize loop info
!
         icloop = -1
         inloop = -1
         ialoop = -1
         itloop = -1
         iriloop = -1
         irfloop = -1
         ictype = 0
         intype = 0
         iatype = 0
         ittype = 0
         iritype = 0
         irftype = 0
         icktype = 0
         bcname = ' '
         lbcname = 1
         baname = ' '
         batag = ' '
         lbaname = 1
         btname = ' '
         lbtname = 1
         biname=bloc_
         mycat=0
         loop_=.false.
         loopnl=0
         nmatch=0
         ksmatch=0
         riname = ' '
         lriname = 0
         rfname = ' '
         lrfname = 0
!
!        Pick up category_keys and list_references
!
         do i = 1,2
210        if(char_(ck(i),name)) then
             if (icktype.ne.0 .and. icktype.ne.i) &
               call warn &
               (' Multiple DDL 1 and 2 related key definitions ')
             icktype = i
             jck = ndict
             call hash_store(locase(name(1:max(1,long_))), &
               dicnam,dicchain, &
               NUMDICT,ndict,dichash,NUMHASH,ick)
             if(ick.eq.0) call err(' CIFdic names > NUMDICT')
             if(ick .eq. jck+1) then
               dictag(ick) = name(1:max(1,long_))
               dictyp(ick) = ' '
               dicxtyp(ick) = ' '
               catkey(ick) = .true.
               alias(ick) = 0
               aroot(ick) = ick
               keychain(ick) = 0
             else
               if(.not.catkey(ick)) then
                 ifind = aroot(ick)
220              catkey(ifind) = .true.
                 ifind = alias(ifind)
                 if (ifind.ne.0) go to 220
               endif
             endif
             if (loop_) go to 210
           endif
         enddo
!
!        Process related items
!
         do i = 1,2
           if(char_(ri(i),name)) then
             if (iritype.ne.0) &
               call warn &
               (' Multiple DDL 1 and 2 related item definitions ')
             iritype = i
             if(loop_) iriloop = loopnl
             riname=locase(name(1:long_))
             lriname=long_
!
!            Seek the matching function, may be in the same loop or not
!
             if(char_(rf(i),name)) then
               if (irftype.ne.0) &
                 call warn &
                 (' Multiple DDL 1 and 2 related item functions ')
               irftype = i
               if (loop_) irfloop = loopnl
               rfname=locase(name(1:long_))
               lrfname=long_
             endif
           endif
         enddo
         loop_ = .false.
         loopnl = 0
!
!        Process categories
!
         do i = 1,3
           if(char_(ct(i),name)) then
             if(ictype.ne.0) &
               call warn(' Multiple DDL 1 and 2 category definitions ')
             ictype = i
             if(loop_) icloop = loopnl
             bcname=locase(name(1:long_))
             lbcname=long_
             nmycat = ndcname+1
             call hash_store(bcname, &
               dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
             if(mycat.eq.0) then
               call err(' Dictionary category names > NUMDICT ')
             endif
             if (mycat.eq.nmycat) then
               ccatkey(mycat) = 0
             endif
!
!            if this is not a loop of categories, we expect a match
!            against the block name, unless we are doing replacements
!
             if(.not.loop_) then
               if(ictype.eq.1) then
                 if(bname(1:lbcname+2).ne. &
                  '_'//bcname(1:lbcname)//'.'  &
                  .and. catchk.eq.'yes' &
                  .and. (rfname(1:7).ne.'replace')) then
                 call warn(' Category id does not match block name')
                 endif
               else
                 if(ictype.eq.2) then
                   if(bcname.ne.'dictionary_definition' .and. &
                      bcname.ne.'category_overview') then
                   if(bname(1:lbcname+2).ne. &
                     '_'//bcname(1:lbcname)//'_') then
                   if(bname(1:lbcname+2).ne. &
                     '_'//bcname(1:lbcname)//' '  &
                  .and. catchk.eq.'yes' &
                  .and. (rfname(1:7).ne.'replace')) then
                   call warn(' Category id does not match block name')
                   endif
                   endif
                   endif
                 endif
               endif
             endif
           endif
           loop_ = .false.
           loopnl = 0
         enddo
!
!        Process names
         do i = 1,2
         if(char_(nt(i),name)) then
           if(intype.ne.0) &
             call warn(' Multiple DDL 1 and 2 name definitions ')
           intype = i
           bxname=locase(name(1:long_))
           if(loop_) inloop = loopnl
         endif
         loop_ = .false.
         loopnl=0
         enddo
         if(intype.eq.0.and.ictype.ne.3.and.(.not.glob_) &
           .and.bname(1:lbname).ne.' ') &
           call warn (' No name defined in block')
         loop_ = .false.
         if(char_(at(1),name)) then
           iatype=1
           baname = locase(name(1:long_))
           batag = name(1:long_)
           lbaname = long_
           if(loop_) ialoop = loopnl
         endif
         loop_ = .false.
         loopnl=0
         if(ictype.ne.3) then
           do i=1,2
             if(char_(tt(i),name)) then
               if(ittype.ne.0) &
                 call warn(' Multiple DDL 1 and 2 type definitions ')
               ittype = i
               btname = locase(name(1:long_))
               if(loop_) itloop = loopnl
             endif
             loop_ = .false.
             loopnl=0
           enddo
         endif
!
!        Now test for consistent combinations
!
         if(inloop.ne.-1) then
           if(icloop.ne.-1.and.icloop.ne.inloop  &
                  .and. catchk.eq.'yes') &
             call warn( &
             ' Categories and names in different loops')
           if(iatype.ne.0.and.ialoop.ne.inloop) then
             if(ialoop.eq.-1) then
               if(bxname.ne.bname) &
                call warn( &
               ' One alias, looped names, linking to first')
             else
               call warn( &
               ' Aliases and names in different loops ' &
               //' only using first alias ')
             endif
           endif
           if(itloop.ne.-1.and.itloop.ne.inloop) &
             call warn( &
             ' Types and names in different loops')
         else
           if(icloop.ne.-1) &
             call warn( &
               ' Multiple categories for one name')
           if(itloop.ne.-1) &
             call warn( &
               ' Multiple types for one name')
         endif
!
!        This is the main loop
!
         if(intype.eq.0) go to 200
250      if(.not.char_(nt(intype),name)) goto 200
         kdict=ndict+1
251      call hash_store(locase(name(1:long_)),dicnam,dicchain, &
           NUMDICT,ndict,dichash,NUMHASH,ifind)
         if(ifind.eq.0) call err(' Cifdic names > NUMDICT')
         if(ifind.eq.kdict) then
           dictag(ifind)=name(1:long_)
           catkey(ifind)=.false.
           aroot(ifind) = ifind
           alias(ifind) = 0
           keychain(ifind) = 0
         endif
         if(ifind.le.idstrt) then
           if(kdup)252,253,254
252        call err(' Duplicate name in dictionary '//dictag(ifind))
253        dicnam(ifind)=char(0)
           goto 251
254        continue       
         endif
         if(dicnam(ifind).eq.bname) nmatch=ifind
         if(dicnam(ifind)(1:lbname).eq.bname) ksmatch=ifind
!dbg     if(dicnam(ifind).ne.bname)
!dbg *   call warn (' Name mismatch: '//dicnam(ifind)//bname)
         if(inloop.ge.0)then
!
!          We are in a loop of names.  If it is the same loop as
!          for categories, we need to extract the matching category
!
           if(inloop.eq.icloop) then
             mycat=0
             if(char_(ct(ictype),name)) then
               bcname=locase(name(1:long_))
               lbcname=long_
               nmycat=ndcname+1
               call hash_store(bcname, &
               dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
               if(mycat.eq.0) then
                 call err(' Dictionary category names > NUMDICT ')
               endif
               if(mycat.eq.nmycat) ccatkey(mycat)=0
             endif
           endif
!
!          If it is the same loop as for types, we need to extract
!          the matching type
!
           if(inloop.eq.itloop) then
             btname=' '
             if(char_(ct(ittype),name)) then
               btname=locase(name(1:long_))
               lbtname=long_
             endif
           endif
!
!          If it is the same loop as for aliases, we need to extract
!          the matching alias
!
           if(inloop.eq.ialoop) then
             baname=' '
             batag=' '
             if(char_(at(1),name)) then
               baname = locase(name(1:long_))
               batag = name(1:long_)
               lbaname = long_
             endif
           endif
         endif
!
!        now we have a name stored in dicnam at location ifind
!        the index of the category in mycat, the type in btname,
!        the alias in baname
!
!        First verify match between the name and category, if
!        we have one, or extract from the block name
!
         if (mycat.eq.0) then
         if (dcindex(ifind).eq.0) then
           if (dicnam(ifind).eq.bloc_) then
             call excat(dicnam(ifind),bcname,lbcname)
!dbg         call warn(' Extracting category name from block name '
!dbg *       //bloc_(1:max(1,lastnb(bloc_))))
             if(bcname(1:1).ne.' ') then
               ictype = 1
               nmycat = ndcname+1
               call hash_store(bcname, &
               dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,mycat)
               if(mycat.eq.0) then
                 call err(' Dictionary category names > NUMDICT ')
               endif
               if (mycat.eq.nmycat) then
                 ccatkey(mycat) = 0
               endif
             else 
               if(catchk.eq.'yes') &
               call warn(' No category defined in block '  &
             //bloc_(1:max(1,lastnb(bloc_)))//' and name ' &
             //dicnam(ifind)(1:max(1,lastnb(dicnam(ifind)))) &
             //' does not match')
             endif
           endif
         endif
         else
         if (bcname(1:lbcname).ne.'dictionary_definition' .and. &
           bcname(1:lbcname).ne.'category_overview') then
           if (dicnam(ifind)(1:lbcname+1).ne.'_'//bcname(1:lbcname) &
              .or.( dicnam(ifind)(lbcname+2:lbcname+2).ne.'_' .and. &
                dicnam(ifind)(lbcname+2:lbcname+2).ne.'.' .and. &
                dicnam(ifind)(lbcname+2:lbcname+2).ne.' ' )) then
                if (catchk.eq.'yes'.and.rfname(1:7).ne.'replace')  &
                call warn(' Item name '// &
                dicnam(ifind)(1:max(1,lastnb(dicnam(ifind))))//' '// &
             ' does not match category name '//bcname(1:lbcname))
           endif
         endif
         endif
!
!        We will need the type in what follows.  cif_mm.dic defines
!        some higher level types.  We map them to primitive types
!
         mapped = btname(1:4)
         do i = 1,12
           if (btname(1:4).eq.map_type(i)) mapped = map_to(i)
         enddo
         if (mapped.ne.'char' .and. &
             mapped.ne.'text' .and. &
             mapped.ne.'    ' .and. &
             mapped.ne.'null' .and. &
             mapped.ne.'numb' ) then
             if (tcheck .eq. 'yes') call warn (' Item type '// &
             btname(1:max(1,lastnb(btname)))//' not recognized')
         endif
!
!        There are two cases to consider, one if the name is new to
!        the dictionary, the other, if it is not
!
         if(ifind.eq.kdict) then
           aroot(ifind)=ifind
           alias(ifind)=0
           dcindex(ifind)=mycat
           dictyp(ifind)=mapped
           dicxtyp(ifind)=btname
         else
           if(dcindex(ifind).ne.mycat) then
             if(dcindex(ifind).eq.0) then
               jfind=ifind
               if (aroot(ifind).ne.0) jfind=ifind
255            continue
               dcindex(jfind)=mycat
               jfind=alias(jfind)
               if(jfind.ne.0) goto 255
             else
               if(mycat.ne.0.and. &
                 (vcheck.eq.'yes'.or.tcheck.eq.'yes') &
                 .and.catchk.eq.'yes')  then
                 if(rfname(1:7).ne.'replace')  &
                 call warn(' Attempt to redefine category for item')
                 endif
             endif
           endif
           if(dictyp(ifind).ne.mapped .or. &
             dicxtyp(ifind).ne.btname) then
             if(dictyp(ifind).eq.' ') then
               jfind=ifind
               if (aroot(ifind).ne.0) jfind=ifind
256            continue
               dictyp(jfind)=mapped
               dicxtyp(jfind)=btname
               jfind=alias(jfind)
               if(jfind.ne.0) go to 256
             else
               if(mapped.ne.' '.and.tcheck.eq.'yes') &
                 call warn(' Attempt to redefine type for item')
             endif
           endif
         endif
!
!        now deal with alias, if any.
!
         if(baname.ne.' ') then
           kadict=ndict+1
           call hash_store(baname(1:lbaname),dicnam,dicchain, &
           NUMDICT,ndict,dichash,NUMHASH,iafind)
           if(iafind.eq.0) call err(' Cifdic names > NUMDICT')
           if(iafind.eq.kadict) then
             dictag(iafind)    =batag
             aroot(iafind)     =aroot(ifind)
             if(aroot(iafind).eq.0) aroot(iafind)=ifind
             catkey(iafind)    =catkey(ifind)
             alias(iafind)     =0
             alias(ifind)      =iafind
             dcindex(iafind)   =dcindex(ifind)
             dictyp(iafind)    =dictyp(ifind)
             dicxtyp(iafind)   =dicxtyp(ifind)
             keychain(iafind)  =0
           else
             if(aroot(iafind).ne.0 .and. &
               aroot(iafind).ne.iafind) then
               if(aroot(iafind).eq.ifind .or. &
                 aroot(iafind).eq.aroot(ifind)) then
                 call warn(' Duplicate definition of same alias')
               else
                 call warn(' Conflicting definition of alias')
               endif
             else
               if((dcindex(iafind).eq.0.or. &
                 dcindex(iafind).eq.dcindex(ifind)).and. &
                 (dictyp(iafind).eq.' '.or. &
                 (dictyp(iafind).eq.dictyp(ifind) .and. &
                  dicxtyp(iafind).eq.dicxtyp(ifind)))) then
                 dcindex(iafind)   =dcindex(ifind)
                 dictyp(iafind)    =dictyp(ifind)
                 dicxtyp(iafind)   =dicxtyp(ifind)
               endif
               aroot(iafind)     =aroot(ifind)
               if(aroot(iafind).eq.0) aroot(iafind)=ifind
               alias(ifind)      =iafind
               if (catkey(iafind)) catkey(ifind) = .true.
               if (catkey(ifind)) catkey(iafind) = .true.
             endif
           endif
         endif
         if(inloop.ge.0) then
           baname = ' '
           batag = ' '
         endif
!
         if(inloop.ge.0.and.loop_) go to 250
         if(nmatch.eq.0) then
         if ((ksmatch.eq.0.or.inloop.lt.0) &
           .and.(rfname(1:7).ne.'replace')) then
         call warn(' No name in the block matches the block name')
         endif
         endif
!
!        check for aliases
!        we execute this loop only in the case of unlooped name
!        with looped alias
!
         if(inloop.lt.0.and.ialoop.ge.0) then
           loop_=.false.
           loopnl=0 
           ganame=baname
260        if(.not.char_(at(iatype),name)) goto 200
           baname=locase(name(1:long_))
           batag=name(1:long_)
           lbaname=long_
           if(baname.eq.ganame) then
             if(loop_) go to 260
             go to 200
           endif
           if(baname.ne.' ') then
             kadict=ndict+1
             call hash_store(baname(1:lbaname),dicnam,dicchain, &
             NUMDICT,ndict,dichash,NUMHASH,iafind)
             if(iafind.eq.0) call err(' CIFdic names > NUMDICT')
             if(iafind.eq.kadict) then
               dictag(iafind)    =batag
               aroot(iafind)     =aroot(ifind)
               if(aroot(iafind).eq.0) aroot(iafind)=ifind
               catkey(iafind)    =catkey(ifind)
               alias(iafind)     =0
               alias(ifind)      =iafind
               dcindex(iafind)   =dcindex(ifind)
               dictyp(iafind)    =dictyp(ifind)
               dicxtyp(iafind)   =dicxtyp(ifind)
               keychain(iafind)  =0
               ifind=iafind
             else
               if(aroot(iafind).ne.0 .and. &
                 aroot(iafind).ne.iafind) then
                 if(aroot(iafind).eq.ifind .or. &
                   aroot(iafind).eq.aroot(ifind)) then
                   call warn(' Duplicate definition of same alias')
                 else
                   call warn(' Conflicting definition of alias')
                 endif
               else
                 if((dcindex(iafind).eq.0.or. &
                 dcindex(iafind).eq.dcindex(ifind)).and. &
                 (dictyp(iafind).eq.' '.or. &
                 (dictyp(iafind).eq.dictyp(ifind) .and. &
                  dicxtyp(iafind).eq.dicxtyp(ifind)))) then
                 dcindex(iafind)   =dcindex(ifind)
                 dictyp(iafind)    =dictyp(ifind)
                 dicxtyp(iafind)   =dicxtyp(ifind)
                 ifind=iafind
                 endif
                 aroot(iafind)     =aroot(ifind)
                 if(aroot(iafind).eq.0) aroot(iafind)=ifind
                 alias(ifind)      =iafind
                 if (catkey(iafind)) catkey(ifind) = .true.
                 if (catkey(ifind)) catkey(iafind) = .true.
               endif
             endif
           endif
           if(loop_) go to 260
         endif
         go to 200
!
400      bloc_=' '
         if (ndcname.ne.0) then
         do ii = idstrt+1,ndict
         if (aroot(ii).eq.0.and.dcindex(ii).eq.0 &
           .and.catchk.eq.'yes') &
           call warn(' No category specified for name '// &
             dicnam(ii)(1:max(1,lastnb(dicnam(ii)))))
         enddo
         endif
         do ii = idstrt+1,ndict
         if (dicxtyp(ii).eq.' ') then
           dicxtyp(ii) = 'null'
           dictyp(ii) = 'null'
           if (tcheck.eq.'yes')  then
             jj = lastnb(dicnam(ii))
             if (jj.gt.0) then
             if (dicnam(ii)(jj:jj).ne.'_') &
             call warn(' No type specified for name '// &
               dicnam(ii)(1:max(1,lastnb(dicnam(ii)))))
             endif
           endif
         endif
         if (catkey(ii)) then
           ifind = aroot(ii)
           mycat = dcindex(ifind)
           if (mycat.ne.0) then
             jj = ccatkey(mycat)
             if (jj.eq.0) then
               ccatkey(mycat) = ifind
             else
410            if (keychain(jj).eq.0) then
                 keychain(jj) = ifind
                 keychain(ifind) = 0
               else
                 if(keychain(jj).ne.ifind) then
                   jj = keychain(jj)
                   goto 410
                 endif
               endif
             endif 
           endif
         endif
         enddo
         if (.not.append_) then
           close(dirdev)
           nrecd=0
         endif
         dictfl='no '
500      continue
         if (append_) then
           nrecd=nrecds
           recend_=recends
           recbeg_=recbegs
         endif
         if(dict_) then
           dicname_=xdicnam
           dicver_ =xdicver
         else
           tcheck = otchk
           vcheck = ovchk
         endif
         if(tcheck.eq.'yes') vcheck='yes'
!dbg     WRITE(6,'(i5,3x,a,2x,a)') (i,dicnam(i),dictyp(i),i=1,ndict)
         return
         end
!
!
!
!
!
! >>>>>> Find position of last non_blank in a string
!
         function lastnb(str)
!
         integer    lastnb
         include   'ciftbx.sys'
         character*(*) str
         integer lenn,ii
         lenn = len(str)
         do 100 ii=lenn,1,-1
         if(str(ii:ii).eq.' ') goto 100
         if(str(ii:ii).ne.tab) goto 120
100      continue
         ii=1
120      lastnb = ii
         return
         end
!
!
!
!
!
! >>>>>> Extract the item.category_id from a save frame name
!
         subroutine excat(sfname,bcname,lbcname)
!
         character*(*) sfname,bcname
         integer lbcname,ii,ic,lastnb,lenn
!
!        Note that this logic works only for item.category_id
!        not for category.id
!
         lenn = lastnb(sfname)
         bcname = ' '
         lbcname = 1
         if (lenn.eq.0.or.sfname(1:1).ne.'_') return
         do ii = 1,lenn-2
         ic = 1+lenn-ii
         if (sfname(ic:ic).eq.'.') then
           bcname = sfname(2:ic-1)
           lbcname = ic-2
           return
         endif
         enddo
         return
         end
!
!
!
!
!
! >>>>>> Open a CIF and copy its contents into a direct access file.
!
         function ocif_(fname)
!
         logical   ocif_
         integer   lastnb
         include  'ciftbx.sys'
         logical   test
         character fname*(*)
         integer   case,i,kp,lp,mp,krpp,mpp
!
         save_=.false.
         glob_=.false.
         jchar=MAXBUF
         lastch=0
         if(line_.gt.MAXBUF) call err(' Input line_ value > MAXBUF')
         if(nrecd.ne.0 .and. (.not.append_)) then
           close(dirdev)
           nrecd=0
           lrecd=0
         endif
!
!        clear the memory resident page buffer
!
         do i = 1,NUMPAGE
         mppoint(i)=0
         enddo
!
         case=ichar('a')-ichar('A')
         tab=char(05)
         if(case.lt.0) goto 100
         tab=char(09)
         bloc_=' '
!
!....... Make sure the CIF is available to open
!
100      file_=fname
         do 120 i=1,MAXBUF
         if(file_(i:i).eq.' ') goto 140
120      continue
140      longf_=i-1
         if (longf_.gt.0) then
           inquire(file=file_(1:longf_),exist=test)
           ocif_=test
           if(.not.ocif_)      goto 200
         else
           file_ = ' '
           longf_ = 1
           ocif_ = .true.
         endif
!
!....... Open up the CIF and a direct access formatted file as scratch
!
         if (file_(1:1).ne.' ') &
         open(unit=cifdev,file=fname,status='OLD',access='SEQUENTIAL', &
                          form='FORMATTED')
         if(nrecd.eq.0) &
         open(unit=dirdev,status='SCRATCH',access='DIRECT', &
                          form='FORMATTED',recl=NUMCPP)
         if(append_ .and. nrecd.ne.0) then
           kp=1
           krpp=NUMCPP/MAXBUF
           lp=(nrecd-1)/krpp+1
           mpp=nrecd-(lp-1)*krpp
           mp=mpp*MAXBUF+1
           mppoint(1)=lp
           if(mp+MAXBUF-1.gt.NUMCPP) then
             mp=1
             lp=lp+1
           else
             read(dirdev,'(a)',rec=lp) pagebuf(kp)
           endif
         else
           kp = 1
           lp = 1
           mp = 1
         endif
!
!....... Copy the CIF to the direct access file
!
160      read(cifdev,'(a,a)',end=180) buffer
         nrecd=nrecd+1
         irecd=nrecd
         if (lastnb(buffer(1:MAXBUF)).gt.line_) &
            call warn(' Input line length exceeds line_')
         pagebuf(kp)(mp:mp+MAXBUF-1) = buffer
         mp = mp+MAXBUF
         if (mp+MAXBUF-1.gt.NUMCPP) then
           write(dirdev,'(a)',rec=lp) pagebuf(kp)
           mppoint(kp)=lp
           lp = lp+1
           kp=kp+1
           if(kp.gt.NUMPAGE) kp=1
           mppoint(kp)=0
           mp=1
         endif
         goto 160
!
180      if(mp.gt.1) then
           pagebuf(kp)(mp:NUMCPP) = ' '
           write(dirdev,'(a)',rec=lp) pagebuf(kp)
           mppoint(kp)=lp
         endif
         lrecd=max(0,recbeg_-1)
         jrecd=max(0,recbeg_-1)
         jrect=-1
         irecd=max(0,recbeg_-1)
         recn_=irecd
         recend_=nrecd
         if (file_(1:1).ne.' ') close(cifdev)
200      return
         end
!
!
!
!
!
! >>>>>> Close off direct access file of the current CIF
!         and reset all data name tables and pointers       
!
         subroutine purge_
!
         include   'ciftbx.sys'
!
         integer i
         if(nrecd.ne.0) close(dirdev)
         do i = 1,NUMPAGE
           mppoint(i)=0
         enddo
         do i = 1,MAXBOOK
           ibkmrk(1,i)=-1
           ibkmrk(2,i)=-1
           ibkmrk(3,i)=-1
           ibkmrk(4,i)=-1
         enddo
         recn_=0
         save_=.false.
         glob_=.false.
         jchar=MAXBUF
         lastch=0
         nrecd=0
         lrecd=0
         irecd=0
         nname=0
         nhash=0
         iname=0
         loopct=0
         loopnl=0
         loop_=.false.
         text_=.false.
         append_=.false.
         recbeg_=0
         recend_=0
         return
         end
!
!
!
!
!
! >>>>>> Store the data names and pointers for the requested data block
!
         function data_(name) 
!
         logical   data_
         logical   wasave
         integer   lastnb
         include  'ciftbx.sys'
         character name*(*),flag*4,temp*(NUMCHAR),ltype*4
         character ctemp*(NUMCHAR)
         character xdname*(NUMCHAR)
         character ydname*(NUMCHAR)
         character locase*(MAXBUF),isbuf*(MAXBUF),lsbuf*(MAXBUF)
         logical   ixcat(NUMDICT)
         integer   ndata,idata,nitem,npakt,i,ii,j,k,kchar,krecd
         integer   jj,icc,idd
         integer   fcatnum,lctemp,isrecd,isjchr,islast
         integer   lsrecd,lsjchr,lslast
         integer   pnname,itpos,ipp,ipj
!
         jchar=MAXBUF
         nname=0
         ndata=0
         nhash=0
         nitem=0
         idata=0
         iname=0
         loopct=0
         loopnl=0
         ltype=' '
         posnam_=0
         posval_=0
         posdec_=0
         posend_=0
         data_=.false.
         wasave=.false.
         loop_=.false.
         text_=.false.
         glob_=.false.
         do ii = 1,MAXBOOK
         ibkmrk(1,ii)=-1
         enddo
         irecd=lrecd
         lrecd=min(nrecd,recend_)
         if(name(1:1).ne.' ') irecd=max(0,recbeg_-1)
         call hash_init(dname,dchain,NUMBLOCK,nname,dhash, &
           NUMHASH)
         call hash_init(cname,cchain,NUMBLOCK,ncname,chash, &
           NUMHASH)
         isrecd=irecd
         isjchr=jchar
         islast=lastch
         lsrecd=isrecd
         lsjchr=isjchr
         lslast=islast
         isbuf=' '
         if(lastch.gt.0)isbuf(1:lastch)=buffer(1:lastch)
         lsbuf=' '
         if(lastch.gt.0)lsbuf(1:lastch)=isbuf(1:lastch)
         xdname=locase(name)
!
!....... Find the requested data block in the file
!
100      lsjchr=isjchr
         call getstr
         isjchr=jchar
         if(irecd.ne.isrecd) then
           lsrecd=isrecd
           lslast=islast
           lsbuf=' '
           if(islast.gt.0)lsbuf(1:islast)=isbuf(1:islast)
           isrecd=irecd
           islast=lastch
           isbuf=' '
           if(lastch.gt.0)isbuf(1:lastch)=buffer(1:lastch)
         endif
         if(type_.eq.'fini')           goto 500
         if(type_.ne.'text')           goto 120
110      call getlin(flag)       
         if(buffer(1:1).ne.';')        goto 110
         jchar=2
         goto 100
120      continue
         if(type_.eq.'save') then
           if(long_.lt.6) then
             if(.not.save_) &
               call err(' Save frame terminator found out of context ')
             wasave=.true.
             save_=.false.
             goto 100
           else
             if(save_) &
               call err(' Prior save frame not terminated ')
             save_=.true.
             if(name.eq.' ')          goto 150
             ydname=locase(strg_(6:long_))
             if(ydname.ne.xdname) goto 100
             goto 150
           endif
         endif
         if(type_.eq.'glob') then
           if(name.ne.' ')            goto 100
           glob_=.true.
           goto 150
         endif
         if(type_.eq.'name'.or.type_.eq.'loop') then
           if(name.ne.' ')            goto 100
           if(.not.wasave) &
             call warn(' Data block header missing ')
           isrecd=lsrecd
           islast=lslast
           isjchr=lsjchr
           isbuf=' '
           if(islast.gt.0)isbuf(1:islast)=lsbuf(1:islast)
           data_=.true.
           bloc_=' '
           itpos=jchar-long_
           if(tabx_) then
           itpos=0
           do ipp=1,jchar-long_
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           endif
           posnam_=itpos
           goto 204           
         endif
         if(type_.ne.'data')          goto 100
         if(name.eq.' ')              goto 150
         ydname=locase(strg_(6:long_))
         if(ydname.ne.xdname)   goto 100
150      data_=.true.
         bloc_=strg_(6:long_)
         itpos=jchar-long_
         if(tabx_) then
         itpos=0
         do ipp=1,jchar-long_
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         endif
         posnam_=itpos
!
!....... Get the next token and identify
!
200      call getstr
!dbg     if(dictfl.eq.'no ')
!dbg *    WRITE(6,*) ltype,type_,loop_,nitem,ndata,idata,iname,nname
!
         if(ltype.ne.'name')                goto 201
         if(type_.eq.'numb')                goto 203
         if(type_.eq.'char')                goto 203
         if(type_.eq.'text')                goto 203
         if(type_.eq.'null')                goto 203
         if(type_.eq.'name'.and.loop_)      goto 204
         call err(' Illegal tag/value construction')
201      if(ltype.ne.'valu')                goto 204
         if(type_.eq.'numb')                goto 202
         if(type_.eq.'char')                goto 202
         if(type_.eq.'text')                goto 202
         if(type_.eq.'null')                goto 202
         goto 204
202      if(nitem.gt.0)                     goto 205
         call err(' Illegal tag/value construction')
203      ltype='valu'
         goto 205
204      ltype=type_
!
205      if(type_.eq.'name')           goto 206
         if(type_.eq.'loop')           goto 210
         if(type_.eq.'data')           goto 210
         if(type_.eq.'save')           goto 210
         if(type_.eq.'glob')           goto 210
         if(type_.ne.'fini')           goto 220
206      if(loop_)                     goto 270
210      if(nitem.eq.0)                goto 215
!
!....... End of loop detected; save pointers
!
         npakt=idata/nitem
         if(npakt*nitem.ne.idata) call err(' Item miscount in loop')
         loopni(loopct)=nitem
         loopnp(loopct)=npakt
         nitem=0
         idata=0
215      if(type_.eq.'name')           goto 270
         if(type_.eq.'data')           goto 300
         if(type_.eq.'save')           goto 300
         if(type_.eq.'glob')           goto 300
         if(type_.eq.'fini')           goto 300
!
!....... Loop_ line detected; incr loop block counter
!
         loop_=.true.
         loopct=loopct+1
         if(loopct.gt.NUMLOOP) call err(' Number of loop_s > NUMLOOP')
         loorec(loopct)=irecd
         loopos(loopct)=jchar-long_
         if(quote_.ne.' ') loopos(loopct)=jchar-long_-1
         itpos=0
         do ipp=1,loopos(loopct)
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         loopox(loopct)=itpos
         goto 200
!
!....... This is the data item; store char position and length
!
220      if(loop_ .and. nitem.eq.0) &
         call err(' Illegal tag/value construction')
         loop_=.false.
!
         i=nname
         if(nitem.gt.0) i=i-nitem+mod(idata,nitem)+1
         if(i.lt.1) call err(' Illegal tag/value construction')
         if(dtype(i).ne.'test')       goto 223
         if(dictfl.eq.'yes')          goto 223
         if(tcheck.eq.'no ')          goto 223
!>>>>    if(long_.eq.1.and.strg_(1:1).eq.'?') goto 223
!>>>>    if(long_.eq.1.and.strg_(1:1).eq.'.') goto 223
         if(type_.eq.'null')          goto 223
         if(type_.eq.'numb')          goto 223
         call warn( ' Numb type violated  '//dname(i))
223      if(nitem.le.0)               goto 224
         idata=idata+1
         if(dtype(i).eq.'null') dtype(i)=type_
         if(dtype(i).eq.'numb' .and. &
           (type_.eq.'char'.or.type_.eq.'text')) dtype(i)='char'
224      if(nname.eq.ndata)           goto 230
         ndata=ndata+1
         if(iloop(ndata).gt.1)        goto 225
         krecd=irecd
         kchar=jchar-long_-1
         if(quote_.ne.' ')kchar=kchar-1
225      continue
         if(dtype(ndata).eq.'    ') dtype(ndata)=type_
         drecd(ndata)=krecd
         dchar(ndata)=kchar
         if(nloop(ndata).gt.0)        goto 230
         nloop(ndata)=0
         iloop(ndata)=long_
!
!....... Skip text lines if present
!
230      if(type_.ne.'text')           goto 200
         if(nloop(ndata).eq.0) dchar(ndata)=0
         if(nloop(ndata).eq.0) iloop(ndata)=long_
250      call getlin(flag)
         if(buffer(1:1).eq.';') then
           jchar=2
           goto 200
         endif
         if(flag.eq.'fini') call err(' Unexpected end of data')
         goto 250
!
!....... This is a data name; store name and loop parameters
!
270      temp=locase(strg_(1:long_))
         k=0
         if(dictfl.ne.'yes' .and. ndict.gt.0) then
           call hash_find(temp, &
             dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,k)
           if(k.ne.0) then
             if(alias_ .and. aroot(k).ne.0) temp=dicnam(aroot(k))
           endif
         endif
         pnname=nname
         call hash_store(temp, &
         dname,dchain,NUMBLOCK,nname,dhash, &
           NUMHASH,j)
         if(j.eq.pnname+1) then
           dtag(j)=strg_(1:long_)
           if(k.ne.0) dtag(j)=dictag(k)
           trecd(j)=irecd
           tchar(j)=jchar-long_
           if(quote_.ne.' ') tchar(j)=jchar-long_-1
           itpos=0
           do ipp=1,tchar(j)
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           xchar(j)=itpos
         endif
         if(j.eq.0) &
           call err(' Number of data names > NUMBLOCK')
         if(k.ne.0)temp=dicnam(k)
         if(j.ne.pnname+1) then
           call warn(' Duplicate data item '// &
           temp(1:max(1,lastnb(temp))))
           goto 200  
         endif
         dtype(nname)=' '
         dxtyp(nname)=' '
         cindex(nname)=0
         ddict(nname)=0
         ctemp='(none)'
         lctemp=6
!
         if(dictfl.eq.'yes' .or. vcheck.eq.'no ') goto 290
         j=k
         if(j.ne.0) then
           ddict(nname)=j
           cindex(nname)=dcindex(j)
           dxtyp(nname)=dicxtyp(j)
           dtype(nname)=dictyp(j)
           if(vcheck.eq.'no ')          goto 280
           if(dictyp(j).eq.'numb') then
             dtype(nname)='test'
           endif
           if(cindex(nname).ne.0) then 
             ctemp=dcname(cindex(nname))
             lctemp=lastnb(ctemp)
             goto 290
           endif   
           goto  280
         endif
         call warn(' Data name '// &
                     temp(1:max(1,lastnb(temp))) &
                     //' not in dictionary!')
280      call excat(temp,ctemp,lctemp)
         if (ctemp.eq.' '.or.'_'//ctemp.eq.temp) then
           ctemp = '(none)'
           lctemp= 6
           if (ndcname.ne.0.and.vcheck.eq.'yes') &
             call warn(' No category defined for ' &
             //temp)
         else
           call hash_find(ctemp, &
             dcname,dcchain,NUMDICT,ndcname,dchash,NUMHASH,j)
           if(j.ne.0) then
             cindex(nname) = j
           else
             ipj=ncname
             call hash_store(ctemp(1:lctemp), &
               cname,cchain,NUMBLOCK,ncname,chash,NUMHASH,j)
             if (j.eq.0) &
               call err(' Number of categories > NUMBLOCK ')
             cindex(nname) = -j
             if (ndcname.gt.0.and.j.eq.ipj+1.and.vcheck.eq.'yes' &
               .and.catchk.eq.'yes') &
               call warn(' Category '// &
               ctemp(1:lctemp)//' first implicitly defined in cif ')
           endif
         endif
!
290      lloop(nname)=0
         nloop(nname)=0
         iloop(nname)=0
         if (nitem.eq.0) fcatnum=cindex(nname)
         if(.not.loop_)               goto 200
         nitem=nitem+1
         if(nitem.gt.NUMITEM) &
           call err(' Items per loop packet > NUMITEM')
         nloop(nname)=loopct
         iloop(nname)=nitem
         if (fcatnum.ne.cindex(nname)) then
           temp = '(none)'
           if (fcatnum.gt.0) temp=dcname(fcatnum)
           if (fcatnum.lt.0) temp=cname(-fcatnum)
           if (ctemp(1:lctemp).ne.temp(1:lastnb(temp)) &
           .and.catchk.eq.'yes') &
           call warn (' Heterogeneous categories in loop '// &
           ctemp(1:lastnb(ctemp))//' vs '// &
           temp(1:lastnb(temp)))
           fcatnum=cindex(nname)
         endif
         goto 200
300      continue
!
!....... Are names checked against dictionary?
!
         if(dictfl.eq.'yes')          goto 500
         if(vcheck.eq.'no '.or.ndict.eq.0) goto 500
         do i=1,nname
           if(dtype(i).eq.'test') dtype(i)='numb'
         enddo
!
!        check for category keys
!
         if(catchk.eq.'yes' .and. ndict.gt.0) then
         do i = 1,ndict
         ixcat(i) = .false.
         enddo
!
!        make a pass marking all used tags and their aliases
!
         do i = 1,nname
           icc=cindex(i)
           idd=ddict(i)
           if(icc.ne.0.and.idd.ne.0) then
             icc = aroot(idd)
310          ixcat(icc) = .true.
             icc = alias(icc)
             if (icc.ne.0) goto 310
           endif
         enddo
!
!        now make a pass making certain the keys are
!        used
!
         do i = 1,nname
           idd=cindex(i)
           if (idd.gt.0) then
             icc=ccatkey(idd)
             if(icc.ne.0) then
             if(aroot(icc).ne.0) icc=aroot(icc)
320          if(icc.ne.0) then
               if(.not.ixcat(icc)) then
                 jj = irecd
                 irecd = drecd(i)
                 call warn(' Category key '// &
                   dictag(icc)(1:lastnb(dictag(icc)))// &
                   ' not given for '// &
                   dcname(idd)(1:lastnb(dcname(idd))))
                 ixcat(icc) = .true.
                 irecd = jj
               endif
               icc = keychain(icc)
               if(icc.ne.0) go to 320
             endif
             endif
           endif      
         enddo
         endif
         
!
!....... End of data block; tidy up loop storage
!
500      lrecd=irecd-1
         if(type_.eq.'save'.and.long_.lt.6) then
           itpos=jchar-long_
           if(tabx_) then
           itpos=0
           do ipp=1,jchar-long_
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           endif
           posval_=itpos
         endif
         irecd=isrecd
         jchar=isjchr
         lastch=islast
         recn_=irecd
         buffer=' '
         if(lastch.gt.0)buffer=isbuf(1:lastch)
         jrecd=irecd
         loop_=.false.
         loopct=0
         if(ndata.ne.nname) call err(' Syntax construction error')
!
!dbg     WRITE(6,'(a)')
!dbg *   ' data name                       type recd char loop leng'
!dbg     WRITE(6,'(a,1x,a,4i5)') (dname(i),dtype(i),drecd(i),dchar(i),
!dbg *              nloop(i),iloop(i),i=1,nname)
!dbg     WRITE(6,'(3i5)') (i,loopni(i),loopnp(i),i=1,loopct)
!
         return
         end
!
!
!
!
!
!
! >>>>>> Get the attributes of data item associated with data name
!
         function test_(temp)
!
         logical   test_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  otestf*3
         character  locase*(MAXBUF)
!
         otestf=testfl
         testfl='yes'
         name=locase(temp)
         test_=.true.   
         if(otestf.eq.'no ')  goto 100
         if(name.eq.nametb)   goto 200
100      call getitm(name)        
200      list_ =loopnl
         if(type_.eq.'null') test_=.false.
         return
         end

!
!
!
!
!
! >>>>>> Set or Reference a bookmark
!
         function bkmrk_(mark)
!
         logical   bkmrk_
         include   'ciftbx.sys'
!
         integer   mark,ii,nitem
         character*4 flag
         bkmrk_=.true.
         if(mark.eq.0) then
           do ii=1,MAXBOOK
             if(ibkmrk(1,ii).lt.0)      goto 100
           enddo
           bkmrk_=.false.
           call warn(' More than MAXBOOK bookmarks requested')
           return
100        mark=ii
           ibkmrk(1,ii)=iname
           ibkmrk(2,ii)=irecd
           ibkmrk(3,ii)=jchar
           if(iname.gt.0) then
             ibkmrk(2,ii) = trecd(iname)
             ibkmrk(3,ii) = tchar(iname)
           endif
           ibkmrk(4,ii)=0
           if(iname.gt.0) then
             if(nloop(iname).ne.0.and. &
               loopnl.eq.nloop(iname).and.loopct.ne.0) then
               nitem=loopni(nloop(iname))
               ibkmrk(2,ii)=looprd(1)
               ibkmrk(3,ii)=max(0,loopch(1)-1)
               ibkmrk(4,ii)=loopct
             endif
           endif
         else
           if(ibkmrk(1,mark).lt.0) then
             bkmrk_=.false.
             return
           endif
           iname=ibkmrk(1,mark)
           irecd=ibkmrk(2,mark)
           loopct=ibkmrk(4,mark)
           loop_=.false.
           text_=.false.
           loopnl=-1
           testfl='no '
           if(iname.gt.0) then
            if(nloop(iname).ne.0.and.loopct.ne.0) then
               nitem=loopni(nloop(iname))
               looprd(nitem+1)=ibkmrk(2,mark)
               loopch(nitem+1)=ibkmrk(3,mark)
               do ii = 1,nitem
                 lloop(ii+iname-iloop(iname))=loopct-1
               enddo
               loopct=loopct-1
               if(lloop(iname).gt.0) then
                 loop_=.true.
                 loopnl=nloop(iname)
               endif
             endif
           endif
           jchar=MAXBUF
           if(irecd.gt.0) then
             irecd=irecd-1
             call getlin(flag)
             jchar=ibkmrk(3,mark)
           endif
           ibkmrk(1,mark)=-1
           mark=0
         endif
         return
         end
!
!
!
!
!
!
! >>>>>> Find the location of the requested item in the CIF
!        The argument "name" may be a data item name, blank
!        for the next such item.  The argument "type" may be
!        blank for unrestricted acceptance of any non-comment
!        string (use cmnt_ to see comments), including loop headers,
!        "name" to accept only the name itself and "valu"
!        to accept only the value, or "head" to position to the
!        head of the CIF.  Except when the "head" is requested,
!        the position is  left after the data item provided.
!
         function find_(name,type,strg)
!
         logical   find_
         include   'ciftbx.sys'
         character  name*(*),type*(*),strg*(*),flag*4
         character  jjbuf*(MAXBUF)
         integer    jjchar,jjrecd,jjlast,jjlrec,jjjrec
!
         find_  = .false.
         strg   = ' '
         long_  = 0
         jjchar = jchar
         jjrecd = lrecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf  = ' '
         if(lastch.gt.0) jjbuf(1:lastch)=buffer(1:lastch)
         if(type.eq.'head') then
           lrecd = min(nrecd,recend_)
           irecd = max(0,recbeg_-1)
           jchar=MAXBUF+1
           call getlin(flag)
           if(flag.eq.'fini')       goto 300
           find_=.true.
           lrecd=max(0,recbeg_-1)
           return
         endif
         if(name.ne.' ') then
           testfl='no '
           call getitm(name)
           if(iname.eq.0) goto 300
           if(type.eq.'valu') then
             list_=loopnl
             strg=strg_(1:long_)
             find_=.true.
             return
           endif
           if(type.eq.'name'.or.loopnl.eq.0) then
             irecd=trecd(iname)-1
             call getlin(flag)
             jchar=tchar(iname)
             posnam_=jchar+1
             call getstr
             strg=strg_(1:long_)
             recn_=irecd
             find_=.true.
             return
           endif
           if(type.eq.' ') then
             irecd=loorec(loopnl)-1
             call getlin(flag)
             jchar=loopos(loopnl)
             call getstr
             posval_=loopos(loopnl)
             if(tabx_) posval_=loopox(loopnl)
             strg=strg_(1:long_)
             recn_=irecd
             find_=.true.
             return
           endif
           call err(' Call to find_ with invalid arguments')
         endif
         if(name.eq.' ') then
200        call getstr
           if(type_.eq.'fini')      goto 300
           if(type.ne.' '.and. &
            (type_.eq.'data'.or.type_.eq.'save'.or. &
            type_.eq.'glob'))   goto 300
           if(type.eq.'name'.and.type_.ne.'name')  goto 200
           if(type.eq.'valu'.and. &
             type_.ne.'numb'.and.type_.ne.'text' &
            .and.type_.ne.'char'.and.type_.ne.'null') goto 200
           find_=.true.
           strg=strg_(1:long_)
           if(type_.eq.'name') then
             posnam_=jchar-long_
           else
             posval_=jchar-long_
             if(quote_.ne.' ') posval_=posval_-1
           endif
           recn_=irecd
           return
         endif
    
!
!        Search failed, restore pointers
!
300      irecd  = jjrecd
         lastch = jjlast
         lrecd  = jjlrec
         jchar  = jjchar
         buffer = ' '
         if(lastch.gt.0)buffer(1:lastch)=jjbuf(1:lastch)
         jrecd  = jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_  = irecd
!           
         return
         end
!
!
!
!
!
!
! >>>>>> Get the next data name in the data block
!
         function name_(temp)
!
         logical    name_
         include   'ciftbx.sys'
         character  temp*(*)
!
         name_=.false.
         temp=' '
         iname=iname+1
         if(iname.gt.nname)  goto 100
         name_=.true.
         temp=dtag(iname)
         if(ddict(iname).ne.0) temp=dictag(ddict(iname))
100      return
         end
!
!
!
!
!
!
! >>>>>> Extract a number data item and its standard deviation
!        This version return single precision numbers
!
         function numb_(temp,numb,sdev)
!
         logical    numb_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  locase*(MAXBUF)
         real       numb,sdev
!
         name=locase(temp)
         if(testfl.eq.'no ')  goto 100
         if(name.eq.nametb)   goto 150
!
100      call getitm(name)
!
150      numb_=.false.
         if(type_.ne.'numb') goto 200
         numb_=.true.
         numb =numbtb
         if(sdevtb.ge.0.0) sdev=sdevtb
!
200      testfl='no '
         return
         end
!
!
!
!
!
!
! >>>>>> Extract a number data item and its standard deviation
!        This version returns double precision numbers
!
         function numd_(temp,numb,sdev)
!
         logical    numd_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  locase*(MAXBUF)
         double precision numb,sdev
!
         name=locase(temp)
         if(testfl.eq.'no ')  goto 100
         if(name.eq.nametb)   goto 150
!
100      call getitm(name)
!
150      numd_=.false.
         if(type_.ne.'numb') goto 200
         numd_=.true.
         numb =numbtb
         if(sdevtb.ge.0.0) sdev=sdevtb
!
200      testfl='no '
         return
         end
!
!
!
!
!
!
! >>>>>> Extract a character data item.
!
         function char_(temp,strg)
!
         logical    char_
         include   'ciftbx.sys'
         character  temp*(*),name*(NUMCHAR)
         character  strg*(*),flag*4
         character  locase*(MAXBUF)
         integer    icpos,itpos,ixpos,ixtpos,ipp,iepos,ispos
!
         name=locase(temp)
         if(testfl.eq.'yes')    goto 100
         if(.not.text_)         goto 120
         if(name.ne.nametb)     goto 120
         char_=.false.
         text_=.false.
         strg=' '
         long_=0
         call getlin(flag)
         if(flag.eq.'fini')    goto 200
         if(buffer(1:1).eq.';') then
           jchar=2
           goto 200
         endif
         quote_=' '
         jchar=lastch+1
         long_=lastch
         strg_(1:long_)=buffer(1:long_)
         goto 150
!
100      if(name.eq.nametb)     goto 150
!
120      call getitm(name)
         if(type_.eq.'null') then
           char_=.false.
           text_=.false.
           strg_=' '
           long_=0
           goto 200
         endif
!
150      char_=.true.
         text_=.false.
         if(tabx_) then
           call detab
           icpos=jchar-long_
           if(quote_.ne.' ') icpos=icpos-1
           iepos=icpos+long_-1
           itpos=0
           do ipp=1,icpos
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
           ispos=itpos
160        ixpos=index(buffer(icpos:iepos),tab)
           ixtpos=itpos+ixpos-1
           if(ixpos.gt.0.and.ixtpos.le.MAXBUF) then
             ixtpos=((ixtpos+7)/8)*8
             icpos=icpos+ixpos
             itpos=ixtpos+1
             if(icpos.le.iepos) goto 160
           else
           strg = &
             bufntb(ispos:min(MAXBUF,itpos+iepos-icpos))
           long_=min(MAXBUF,itpos+iepos-icpos)-ispos+1
           if(ispos.eq.1.and.strg(1:1).eq.';') &
             strg(1:1) = ' '
           endif
         else
           strg=' '
           if(long_.gt.0) then
             strg=strg_(1:long_)
           endif
         endif
         if(type_.eq.'char')   goto 200
         char_=.false.
         if(type_.ne.'text')   goto 200
         char_=.true.
         call getlin(flag)
         jchar=MAXBUF+1
         if(flag.eq.'fini')    goto 200
         if(buffer(1:1).eq.';')then
           jchar=2
           goto 200
         endif
         irecd=irecd-1
         text_=.true. 
!
200      testfl='no '
         return
         end
!
!
!
!
!
!
! >>>>>> Extract a comment field.               
!
         function cmnt_(strg)          
!
         logical   cmnt_
         integer   lastnb
         include  'ciftbx.sys'
         character strg*(*),flag*4,c*1, &
           jjbuf*(MAXBUF)
         integer   jjchar,jjrecd,jjlast,jjlrec,jjjrec
         integer   ipp,itpos,ixpos
!
         jjchar = jchar
         jjrecd = irecd
         jjlast = lastch
         jjlrec = lrecd
         jjjrec = jrecd
         jjbuf=' '
         if(lastch.gt.0)jjbuf(1:lastch)=buffer(1:lastch)
         lrecd = nrecd
         if(bloc_.eq.' ') then
           if(irecd.eq.0) jchar=MAXBUF
         endif
         strg=' '
         long_=0
         cmnt_=.false.
         goto 105
100      jchar=jchar+1
105      if(jchar.le.lastch)     goto 140
!
!....... Read a new line
!
110      call getlin(flag)
         if(flag.eq.'fini') then
           strg='fini'
           jchar=MAXBUF+1
           long_=4
           cmnt_=.false.
           return
         endif
         jchar=1
         strg=char(0)
         long_=1
         posnam_=0
         goto 220
140      if(lastch.eq.1.and.buffer(1:1).eq.' ') go to 200
!
!....... Process this character in the line
!
150      c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab.and.(.not.tabx_)) goto 190
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 200
         goto 300
!
!        For a tab, when not expanding to blanks, accept
!        that single character as a comment
!
190      long_=1
         strg=tab
         posnam_=jchar
         jchar=jchar+1
         goto 220
!
!....... Accept the remainder of the line as a comment 
!
200      long_=lastch-jchar
         itpos=jchar
         if(tabx_) then
           itpos=0
           do ipp=1,jchar
             itpos=itpos+1
             if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
           enddo
         endif
210      posnam_=itpos
         if(long_.gt.0) then
           if(tabx_) then
             call detab
             ixpos= lastnb(bufntb)
             strg = bufntb(itpos+1:ixpos)
           else
             strg = buffer(jchar+1:lastch)
           endif
         endif
         if(long_.le.0) then
           strg=' '
           long_=1
         endif
         jchar=MAXBUF+1
220      lrecd=jjlrec
         cmnt_=.true.
         return
!
!....... Found a non-comment field, restore pointers
!
300      irecd = jjrecd
         lastch = jjlast
         lrecd = jjlrec
         jchar = jjchar
         buffer=' '
         if(lastch.gt.0)buffer(1:lastch)=jjbuf(1:lastch)
         jrecd=jjjrec
         if(jrecd.ne.irecd) jrecd=-1
         recn_=irecd
         return
         end
!
!
!
!
!
! >>>>> Convert name string to lower case
!        
         function locase(name)
!
         include     'ciftbx.sys'
         character    locase*(MAXBUF)
         character    temp*(MAXBUF),name*(*)
         character    low*26,cap*26,c*1
         integer i,j
         data  cap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
         data  low /'abcdefghijklmnopqrstuvwxyz'/
!
         temp=name
         do 100 i=1,MAXBUF
         c=temp(i:i)
         if(c.eq.' ') goto 200
         if(c.eq.tab) goto 200
         j=index(cap,c)
         if(j.ne.0) temp(i:i)=low(j:j)        
100      continue
200      locase=temp
         return
         end
!
!
!
!
!
! >>>>>> Get the data item associated with the tag.
!
         subroutine getitm(name)
!
         include   'ciftbx.sys'
         SAVE
         character name*(*)
         character flag*4
         integer   iitem,nitem,npakt
         integer   kchar,loopi,i,j,itpos,ipp
!
!....... Find requested dataname in hash list
!
         nametb=name
         posnam_=0
         posval_=0
         posdec_=0
         posend_=0
         quote_=' '
         if(name(1:1).eq.'_')       goto 100
         type_='null'
         dictype_='null'
         diccat_='(none)'
         dicname_=name
         tagname_=' '
         strg_=' '
         long_=1
         goto 1000
100      call hash_find(nametb, &
           dname,dchain,NUMBLOCK,nname,dhash,NUMHASH, &
           iname)
         if(iname.gt.0)             goto 180
         if(dictfl.ne.'yes') then
         call hash_find(nametb, &
           dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,j)
         if(j.ne.0) then
           dictype_=dicxtyp(j)
           if(dcindex(j).ne.0) diccat_=dcname(dcindex(j))
           dicname_=nametb
           if(aroot(j).ne.0) then
             dicname_=dictag(aroot(j))
             call hash_find(dicnam(aroot(j)), &
               dname,dchain,NUMBLOCK,nname,dhash,NUMHASH, &
               iname)
             if(iname.gt.0)      goto 180
           endif
           type_='null'
           tagname_=' '
           strg_=' '
           long_=1
           go to 1000
         endif
         endif
160      continue
         type_='null'
         dictype_='null'
         diccat_='(none)'
         dicname_=name
         long_=1
         goto 1000
!
!
180      tagname_=dtag(iname)
         if(ddict(iname).ne.0) tagname_=dictag(ddict(iname))
         posnam_=tchar(iname)
         if(tabx_)posnam_=xchar(iname)
         if(nloop(iname).le.0)      goto 500
!
!....... Process loop packet if first item request
!
         if(nloop(iname).ne.loopnl) goto 200
         if(lloop(iname).lt.loopct) goto 300
         if(loop_)                  goto 230
200      loop_=.true.
         loopct=0
         loopnl=nloop(iname)
         nitem=loopni(loopnl)
         npakt=loopnp(loopnl)
         irecd=drecd(iname)-1
         call getlin(flag)
         jchar=max(0,dchar(iname)-1)
!dbg     if(jchar.lt.0) write(6,'(7H dchar ,i5)') jchar
         do 220 i=1,nitem
220      lloop(i+iname-iloop(iname))=0
         goto 240
!
!....... Read a packet of loop items
!
230      nitem=loopni(loopnl)
         npakt=loopnp(loopnl)
         irecd=looprd(nitem+1)-1
         call getlin(flag)
         jchar=loopch(nitem+1)
!dbg     if(jchar.lt.0) write(6,'(7H loopch,i5)') jchar
240      iitem=0
250      iitem=iitem+1
         if(iitem.le.nitem)     goto 255
         loopch(iitem)=jchar
         looprd(iitem)=irecd
         goto 270
255      call getstr
         loopch(iitem)=jchar-long_
         if(quote_.ne.' ')loopch(iitem)=jchar-long_-1
         loopln(iitem)=long_
         looprd(iitem)=irecd
         if(buffer(1:1).ne.';'.or.loopch(iitem).ne.1) &
                                goto 250
260      call getlin(flag)
         if(flag.eq.'fini') call err(' Unexpected end of data')
         if(buffer(1:1).ne.';') goto 260
         jchar=2
         goto 250
270      loopct=loopct+1
         if(loopct.lt.npakt)    goto 300
         loop_=.false.
!
!....... Point to the loop data item
!
300      lloop(iname)=lloop(iname)+1
         loopi=iloop(iname)
         irecd=looprd(loopi)-1
         call getlin(flag)
         long_=loopln(loopi)
         kchar=loopch(loopi)
         goto 550
!
!....... Point to the non-loop data item
!
500      irecd=drecd(iname)-1
         call getlin(flag)
         kchar=dchar(iname)+1
         long_=iloop(iname)
         loop_=.false.
         loopct=0
         loopnl=0
!
!....... Place data item into variable string and make number
!
550      type_=dtype(iname)
         dictype_=dxtyp(iname)
         diccat_='(none)'
         if(cindex(iname).gt.0) diccat_=dcname(cindex(iname))
         if(cindex(iname).lt.0) diccat_=cname(-cindex(iname))
         if(diccat_.eq.' ') diccat_='(none)'
         dicname_=dtag(iname)
         if(ddict(iname).ne.0) then
           if (aroot(ddict(iname)).ne.0) then
             dicname_=dictag(aroot(ddict(iname)))
           endif
         endif
         strg_=' '
         if(long_.gt.0) then
           strg_(1:long_)=buffer(kchar:kchar+long_-1)
         endif
         itpos=kchar
         if(tabx_) then
         itpos=0
         do ipp=1,kchar
           itpos=itpos+1
           if(buffer(ipp:ipp).eq.tab) itpos=((itpos+7)/8)*8
         enddo
         endif
         posval_=itpos
         posend_=itpos+long_-1
         jchar=kchar+long_
         if(jchar.le.MAXBUF) then
           if(buffer(jchar:jchar).ne.' ' .and. &
             buffer(jchar:jchar).ne.tab) jchar=jchar+1
         endif
         quote_=' '
         if(kchar.gt.1) then
           if(buffer(kchar-1:kchar-1).ne.' ' .and. &
              buffer(kchar-1:kchar-1).ne.tab) then
             quote_=buffer(kchar-1:kchar-1)
           endif
         endif
         if(type_.eq.'char' .and. kchar.eq.1 .and. &
           buffer(1:1).eq.';') type_='text'
         if(type_.eq.'text') then
           if(buffer(1:1).eq.';') then
             strg_(1:1)=' '
           else
             type_='char'
           endif
         endif
         if(type_.eq.'numb') then
           call ctonum
           if(posdec_.gt.0) posdec_=posval_+posdec_-1
         endif
         if(type_.eq.'char' .and. strg_.eq.' '.and.nblank_)  &
           type_='null'
         if(quote_.ne.' ') goto 1000
         if(long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
         if(long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
!
1000     return
         end
!
!
!
!
!
!
!
! >>>>>> Read the next string from the file
!
!
         subroutine getstr
!
!        On entry, jchar is set to one less than the next character
!        to be read, on the line given by irecd, which is assumed
!        to have been loaded into buffer, with lastch set to the
!        position of the last character
!
         include   'ciftbx.sys'
         integer   i,j,jj(11),im
         logical   quoted
         character c*1,num*21,flag*4
         data num/'0123456789+-.()EDQedq'/
!
         quoted=.false.
         quote_=' '
         if(irecd.gt.0.and. &
           jchar.le.1.and.lastch.gt.0) then
           jchar=1
           goto 140
         end if
100      jchar=jchar+1
         if(jchar.le.lastch)     goto 150
!
!....... Read a new line
!
110      call getlin(flag)
         type_='fini'
         dictype_=type_
         diccat_='(none)'
         dicname_=' '
!dbg     write(6,'(/5i5,a)') 
!dbg *              irecd,jrecd,lrecd,nrecd,lastch, buffer(1:lastch)
         if(flag.eq.'fini')  goto 500
!
!....... Test if the new line is the start of a text sequence
!
140      if(buffer(1:1).ne.';') goto 150
         type_='text'
         jchar=lastch+1
         long_=lastch
         strg_(1:long_)=buffer(1:long_)
         strg_(1:1)=' '
         goto 500
!
!....... Process this character in the line
!
150      c=buffer(jchar:jchar)
         if(c.eq.' ')       goto 100
         if(c.eq.tab)       goto 100
         if(c.eq.'#')       goto 110
         if(c.eq.'''')      goto 300
         if(c.eq.'"')       goto 300
         if(c.ne.'_')       goto 200
         type_='name'
         goto 210
!
!....... Span blank delimited token; test if a number or a character 
!
200      type_='numb'
         im=0
         do 205 i=1,11
205      jj(i)=0
210      do 250 i=jchar,lastch
         if(buffer(i:i).eq.' ')       goto 400
         if(buffer(i:i).eq.tab)       goto 400
         if(type_.ne.'numb')          goto 250
         j=index(num,buffer(i:i))
         if(j.eq.0)                 type_='char'
         if(j.le.10) then
           im=im+1
           goto 250
         endif
         if(j.gt.13.and.im.eq.0) type_='char'
         jj(j-10)=jj(j-10)+1
250      continue
         i=lastch+1
         if(type_.ne.'numb') goto 400
         do 270 j=1,5
         if((jj(j).gt.1.and.j.gt.2) .or. &
           jj(j).gt.2)             type_='char'
270      continue
         goto 400
!
!....... Span quote delimited token; assume character
!
300      type_='char'
         quoted=.true.
         jchar=jchar+1
         do 320 i=jchar,lastch
         if(buffer(i:i).ne.c)             goto 320
         if(i+1.ge.lastch)                goto 400
         if(buffer(i+1:i+1).eq.' ')       goto 400
         if(buffer(i+1:i+1).eq.tab)       goto 400
320      continue
!dbg     write(6,'(a,4i5,a)') 
!dbg *       '**** ',irecd,lastch,i,jchar,buffer(jchar:i)       
         call warn(' Quoted string not closed')
!
!....... Store the string for the getter
!
400      long_=0
         strg_=' '
         if(i.gt.jchar) then
           long_=i-jchar
           strg_(1:long_)=buffer(jchar:i-1)
         endif
         jchar=i
         quote_=' '
         if(quoted) then
           quote_=buffer(jchar:jchar)
           jchar =jchar+1
         endif
!dbg     write(6,'(5x,8i5,5x,a)') 
!dbg *   irecd,jrecd,lrecd,nrecd,lastch,i,jchar,long_,strg_(1:long_)
         if(type_.ne.'char'.or.quoted) goto 500
         if(strg_(1:5).eq.'data_') type_='data'
         if(strg_(1:5).eq.'loop_') type_='loop'
         if(long_.eq.1.and.strg_(1:1).eq.'?') type_='null'
         if(long_.eq.1.and.strg_(1:1).eq.'.') type_='null'
         if(strg_(1:5).eq.'save_') type_='save'
         if(long_.eq.7.and. strg_(1:7).eq.'global_') type_='glob'
!
500      return
         end
!
!
!
!
!
!
! >>>>>> Convert a character string into a number and its esd
!
!                                          Q
!                                          D+
!                                          E-
!                                +         +
!           number string        -xxxx.xxxx-xxx(x)
!           component count CCNT 11111222223333444
!           (with at least 1 digit in the mantissa)
!
         subroutine ctonum
!
         integer   lastnb
         include  'ciftbx.sys'
         character test*22,c*1
         integer*4 m,nchar
         integer*4 ccnt,expn,msin,esin,ndec,ids,nmd
         double precision numb,sdev,ntemp,mant
         data test /'0123456789+.-()EDQedq '/
!
         numbtb=0.D0
         sdevtb=-1.D0
         numb=1.D0
         sdev=0.D0
         ccnt=0
         mant=0.D0
         expn=0.
         msin=+1
         esin=+1
         ndec=0
         ids=0
         nmd=0
         type_='char'
         posdec_=0
         esddig_=0
         if(long_.eq.1.and. &
           index('0123456789',strg_(1:1)).eq.0) goto 500
         lzero_=.false.
         decp_=.false.
!
!....... Loop over the string and identify components
!
!        The scan works in phases
!          ccnt = 0   processing looking for first digit
!          ccnt = 1   processing before decimal point
!          ccnt = 2   processing after decimal point
!          ccnt = 3   processing exponent
!          ccnt = 4   processing standard deviation
!
         do 400 nchar=1,long_
!
         c=strg_(nchar:nchar)
         m=index(test,c)
         if(m.eq.0)     goto 500
         if(m.gt.10)    goto 300
!
!....... Process the digits
!
         if(ccnt.eq.0)  ccnt=1
         if(ccnt.eq.2)  ndec=ndec+1
         if(ccnt.gt.2)  goto 220
         ntemp=m-1
         if (ndec.eq.0) then
           mant=mant*10.D0+ntemp
         else
           mant=mant+ntemp/10.D0**(ndec)
         endif
         nmd=nmd+1
         if(ccnt.eq.1.and.mant.ne.0.D0) ids=ids+1
         goto 400
220      if(ccnt.gt.3)  goto 240
         expn=expn*10+m-1
         goto 400
240      esddig_=esddig_+1
         ntemp=m-1
         sdev=sdev*10.D0+ntemp
         sdevtb=1.D0
         goto 400
!
!....... Process the characters    . + - ( ) E D Q
!
300      if(c.ne.'.')  goto 320
         decp_=.true.
         if(nchar.gt.1.and.mant.eq.0.d0) then
           if(strg_(nchar-1:nchar-1).eq.'0') lzero_=.true.
         endif
         if(ccnt.gt.1) goto 500
         posdec_=nchar
         ccnt=2
         goto 400
!
320      if(nmd.eq.0.and.m.gt.13) goto 500
         if(c.ne.'(')  goto 340
         if(posdec_.eq.0) posdec_=nchar
         ccnt=4
         goto 400
!
340      if(posdec_.eq.0.and.ccnt.gt.0) posdec_=nchar
         if(c.eq.' ')  goto 400
         if(m.gt.13) m = 11
         if(ccnt.eq.3) goto 500
         if(ccnt.gt.0) goto 360
         ccnt=1
         msin=12-m
         goto 400
360      ccnt=3
         esin=12-m
!
400      continue
!
         if(posdec_.eq.0) posdec_=lastnb(strg_(1:long_))+1
!
!....... String parsed; construct the numbers
!
         expn=expn*esin
         if(expn+ids.gt.-minexp) then
           call warn(' Exponent overflow in numeric input')
           expn=-minexp-ids
         endif
         if(expn.lt.minexp) then
           call warn(' Exponent underflow in numeric input')
           expn=minexp
         endif
         if(expn-ndec.lt.0) numb=1./10.D0**abs(expn-ndec)
         if(expn-ndec.gt.0) numb=10.D0**(expn-ndec)
         if(sdevtb.gt.0.0) sdevtb=numb*sdev
         numb=1.D0
         if(expn.lt.0) numb=1./10.D0**abs(expn)
         if(expn.gt.0) numb=10.D0**(expn)
         ntemp=msin
         numbtb=numb*mant*ntemp
         type_='numb'
!
500      return
         end
!
!
!
!
!
!
! >>>>>> Read a new line from the direct access file
!
         subroutine getlin(flag)
!
         integer   lastnb
         include  'ciftbx.sys'
         character flag*4
         integer   krpp,kpp,lpp,mpp,npp,ir
!
         irecd=irecd+1
         jchar=1
         if(irecd.eq.jrecd.and. &
           irecd.gt.recbeg_.and. &
           irecd.le.recend_)  goto 200
         if(irecd.le.min(lrecd,recend_))  goto 100
         irecd=min(lrecd,recend_)+1
         buffer=' '
         lastch=0
         jchar=MAXBUF+1
         jrecd=-1
         flag='fini'
         goto 200
100      continue
         lpp=-1
         mpp=-1
         npp=kpp
         krpp=NUMCPP/MAXBUF
         kpp=(irecd-1)/krpp+1
         do ir = 1,NUMPAGE
           if(mppoint(ir).eq.kpp) then
             lpp = ir
             goto 120
           endif
           if(mppoint(ir).eq.0) then
             lpp=ir
           else
             if(iabs(mppoint(ir)-kpp) &
               .gt.iabs(npp-kpp)) then
               mpp=ir
               npp=mppoint(ir)
             endif
           endif
         enddo
!
!        failed to find page as resident
!        remove a target page
!
         if(lpp.eq.-1)lpp=mpp
         if(lpp.eq.-1)lpp=1
         mppoint(lpp)=kpp
         read(dirdev,'(a)',rec=kpp) pagebuf(lpp)
120      mpp=irecd-1-(kpp-1)*krpp
         npp=mpp*MAXBUF+1
         buffer=pagebuf(lpp)(npp:npp+MAXBUF-1)
         recn_=irecd
         lastch=max(1,lastnb(buffer))
         jrecd=irecd
         flag=' '
200      return
         end
!
!
!
!
!
!
! >>>>>> Detab buffer into bufntb
!
         subroutine detab
!
         include   'ciftbx.sys'
         integer   icpos,itpos,ixpos,ixtpos
         if(jrecd.eq.jrect) return
         icpos=1
         itpos=1
         bufntb=' '
         if(lastch.gt.0) then
100      ixpos=index(buffer(icpos:lastch),tab)
         ixtpos=ixpos+itpos-1
         if(ixpos.gt.0.and.ixtpos.le.MAXBUF) then
           ixtpos=((ixtpos+7)/8)*8
           if(ixpos.gt.1) then
           bufntb(itpos:ixtpos)= &
             buffer(icpos:ixpos+icpos-2)
           else
           bufntb(itpos:ixtpos)=' '
           endif
           itpos=ixtpos+1
           icpos=ixpos+icpos
           goto 100
         else
           bufntb(itpos:max(MAXBUF,itpos+lastch-icpos))= &
             buffer(icpos:lastch)
         endif
         endif
         jrect=jrecd
         return
         end
!
!
!
!
!
!
! >>>>>> Write error message and exit.
!
         subroutine err(mess)
         character*(*) mess
         call cifmsg('error',mess)
         stop
         end
!
!
!
!
!
!
! >>>>>> Write warning message and continue.
!
         subroutine warn(mess)
         character*(*) mess
         call cifmsg('warning',mess)
         return
         end
!
!
!
!
!
!
! >>>>>> Write a message to the error device
!
         subroutine cifmsg(flag,mess)
!
         integer    lastnb
         include   'ciftbx.sys'
         character*(*) flag
         character*(*) mess
         character*(MAXBUF)  tline
         character*5   btype
         integer       ll,ls,ltry,ii,i
!
         btype = 'data_'
         if(save_) btype = 'save_'
         if(.not.glob_) then
         tline= ' ciftbx '//flag//': ' &
         //file_(1:longf_)//' '//btype &
         //bloc_(1:max(1,lastnb(bloc_)))//' line:'
         else
         tline= ' ciftbx '//flag//': ' &
         //file_(1:longf_)//' global_'//' line:'
         endif
         ll = max(1,lastnb(tline))
         write(errdev,'(a,i7)')tline(1:ll),irecd
         ll=len(mess)
         ls=1
100      if(ll-ls.le.79) then
           write(errdev,'(1X,a)') mess(ls:ll)
           return
         else
           ltry = min(ll,ls+79)
           do ii = ls+1,ltry
           i = ltry-ii+ls+1
           if(mess(i:i).eq.' ') then
             write(errdev,'(1X,a)') mess(ls:i-1)
             ls=i+1
             if(ls.le.ll) go to 100
             return
           endif
           enddo
           write(errdev,'(1X,a)') mess(ls:ltry)
           ls=ltry+1
           if(ls.le.ll) go to 100
           return
         endif  
         end
!
!
!
!
! >>>>>> Create a named file.
!
         function pfile_(fname)
!
         logical   pfile_
         include   'ciftbx.sys'
         logical   test
         integer   i
         character fname*(*)
!
!....... Test if a file by this name is already open.
!
         if(pfilef.eq.'yes') call close_
         pfilef='no '
         file_=fname
         do 120 i=1,MAXBUF
         if(file_(i:i).eq.' ') goto 140
120      continue
140      if (i.gt.1) then
           inquire(file=file_(1:i-1),exist=test)
           pfile_=.false.
           longf_ = i-1
           if(test)            goto 200
         else
           file_ = ' '
           pfile_ = .true.
           longf_ = 1
         endif
!
!....... Open up a new CIF
!
         if (file_(1:1) .ne. ' ')  then 
         open(unit=outdev,file=fname,status='NEW',access='SEQUENTIAL', &
                          form='FORMATTED')
         precn_=0
         endif
         pfile_=.true.  
         pfilef='yes'
         nbloc=0
         pchar=1+lprefx
         pcharl=0
         obuf=prefx
         obuf(pchar:MAXBUF)=' '
200      return
         end
!
!
!
!
!
! >>>>>> Store a data block command in the CIF
!        Call with blank name to close current block only
!
         function pdata_(name) 
!
         logical   pdata_
         SAVE
         include  'ciftbx.sys'
         character name*(*),temp*(MAXBUF)
         character dbloc(100)*(NUMCHAR)
         integer   i
!
         pdata_=.true.
         if(ploopn.ne.0)     call eoloop
         if(ptextf.eq.'yes') call eotext
         if(psaveo) then
           pchar=-1
           if(pposval_.ne.0) then
             pchar=lprefx+1
             call putstr(' ')
             pchar=lprefx+pposval_
             pposval_=0
           endif
           call putstr('save_')
           psaveo=.false.
         endif
         if(globo_) then
           pchar=-1
           temp='global_'
           psaveo=.false.
           goto 135
         endif
!
!....... Check for duplicate data name
!
         temp=name
         if(temp.eq.' ')        goto 200
         if(saveo_)             goto 130
         pdata_=.false.
         do 110 i=1,nbloc
         if(temp.eq.dbloc(i))   goto 130
110      continue
         pdata_ = .true.
         goto 125
!
!....... Save block name and put data_ statement
!
125      nbloc=nbloc+1
         if(nbloc.le.100) dbloc(nbloc)=temp
130      pchar=-1
         temp='data_'//name
         if(saveo_) temp='save_'//name
         if(globo_) temp='global_'
         psaveo=saveo_
135      if(pposnam_.gt.0) then
           pchar=lprefx+1
           call putstr(' ')
           pchar=lprefx+pposnam_
           pposnam_=0
         endif
         call putstr(temp)         
         pchar=lprefx
!
200      return
         end
!
!
!
!
!
!
! >>>>>> Put a number into the CIF, perhaps with an esd appended
!
         function pnumb_(name,numb,sdev)
!
         logical    pnumb_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
         real       numb,sdev
         double precision dnumb,dsdev,dprec
!
         pnumb_=.true.
         flag  =.true.
         tflag =.true.
         temp=name
         if(ptextf.eq.'yes') call eotext
!
         if(name(1:1).eq.' ')   goto 120
         if(vcheck.eq.'no ')    goto 100
         call dcheck(temp,'numb',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0) &
             temp=dictag(aroot(xdchk))
         endif
         pnumb_=flag
100      if(ploopn.ne.0)        call eoloop
         pchar=-1
         if(pposnam_.ne.0)pchar=pposnam_+lprefx
         call putstr(temp)
!
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         dprec=decprc
         dnumb=numb
         dsdev=sdev
         call putnum(dnumb,dsdev,dprec)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not in dictionary')
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not correct type')
         endif
!
150      pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         pesddig_=0
         return
         end
!
!
!
!
!
!
! >>>>>> Put a double precision number into the CIF, perhaps 
!        with an esd appended
!
         function pnumd_(name,numb,sdev)
!
         logical    pnumd_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
         double precision numb,sdev
!
         pnumd_=.true.
         flag  =.true.
         tflag =.true.
         temp=name
         if(ptextf.eq.'yes') call eotext
!
         if(name(1:1).eq.' ')   goto 120
         if(vcheck.eq.'no ')    goto 100
         call dcheck(temp,'numb',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0) &
             temp=dictag(aroot(xdchk))
         endif
         pnumd_=flag
100      if(ploopn.ne.0)        call eoloop
         pchar=-1
         if(pposnam_.ne.0)pchar=pposnam_+lprefx
         call putstr(temp)
!
120      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         call putnum(numb,sdev,dpprc)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not in dictionary')
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not correct type')
         endif
!
150      pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         pesddig_=0
         return
         end
!
!
!
!
!
!
! >>>>>> Put a character string into the CIF.
!
         function pchar_(name,string)      
!
         logical    pchar_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR),string*(*)
         character  line*(MAXBUF),strg*(MAXBUF)
         integer    i,j
!
         pchar_=.true.
         flag  =.true.
         tflag =.true.
         temp  =name
         if(ptextf.eq.'yes') call eotext
!
         if(name(1:1).eq.' ')   goto 110
         if(vcheck.eq.'no ')    goto 100
         call dcheck(temp,'char',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0) &
             temp=dictag(aroot(xdchk))
         endif
         pchar_=flag
100      if(ploopn.ne.0)        call eoloop
         pchar=-1
         if(pposnam_.gt.0) pchar=posnam_+lprefx
         call putstr(temp)
!
110      if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         line=string
         do 120 i=MAXBUF,2,-1
         if(line(i:i).ne.' ') goto 130
120      continue
130      if(pposval_.ne.0.and.pposend_.ge.pposval_) &
            i=max(i,pposend_-pposval_+1)
         if(pquote_.ne.' ')   goto 150
         do 140 j=i,1,-1
         if(line(j:j).eq.' ') goto 150
140      continue
         if((line(1:1).eq.'_'  &
           .or. line(i:i).eq.'_' &
           .or. line(1:1).eq.'''' &
           .or. line(1:1).eq.'"' &
           .or. line(1:1).eq.';') &
           .and.line(1:i).ne.'''.''' &
           .and.line(1:i).ne.'''?''' &
           .and.line(1:i).ne.'"."' &
           .and.line(1:i).ne.'"?"') goto 150
         strg=line(1:i)
         goto 200
150      if(pquote_.eq.';')       goto 190
         if(line(1:i).eq.' '.and.nblanko_) then
           strg = '.'
           i = 1
           if(pposval_.ne.0) then
             pchar=pposval_+lprefx
           endif
           call putstr(strg(1:i))   
           go to 210
         endif
         if(pquote_.eq.'"')       goto 170
         do 160 j=1,i-1
         if(line(j:j).eq.''''.and. &
           (line(j+1:j+1).eq.' '.or.line(j+1:j+1).eq.tab)) &
           goto 170
160      continue
165      strg=''''//line(1:i)//''''
         i=i+2
         pquote_=''''
         goto 200
170      do 180 j=1,i-1
         if(line(j:j).eq.'"'.and. &
           (line(j+1:j+1).eq.' '.or.line(j+1:j+1).eq.tab)) &
           goto 190
180      continue
185      strg='"'//line(1:i)//'"'
         i=i+2
         pquote_='"'
         goto 200
190      pchar=-1
         strg='; '//line(1:i)
         i=i+2
         ptextf='yes'
         call putstr(strg(1:i))
         pchar=-1
         ptextf='no '
         call putstr(';')
         pchar=lprefx
         call putstr(' ')
         call warn(' Converted pchar_ output to text for: ' &
           //strg(3:i))
         goto 210
!
200      if(pposval_.ne.0) then
           pchar=pposval_+lprefx
           if(pquote_.ne.' ') pchar=pchar-1
         endif
         call putstr(strg(1:i))   
210      if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not in dictionary')
         endif
         if((.not.tflag).and.line(1:i).ne.'.'.and. &
           line(1:i).ne.'?'.and.pquote_.eq.' ') then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not correct type')
         endif
250      pposval_=0
         pposdec_=0
         pposnam_=0
         pposend_=0
         pquote_=' '
         return
         end
!
!
!
!
!
! >>>>>> Put a comment in the output CIF 
!
         function pcmnt_(string)     
!
         logical    pcmnt_
         include   'ciftbx.sys'
         character  string*(*), temp*(MAXBUF)
!
         if(ptextf.eq.'yes') call eotext
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         if(string.eq.' '.or. &
           (string.eq.char(0)) .or. &
           (string.eq.tab.and.(.not.ptabx_))) then
           if(string.eq.' ') pchar=-1
           call putstr(string)
           if(string.eq.' ') call putstr(char(0))
         else
           temp='#'//string
           call putstr(temp)
           call putstr(char(0))
         endif
         pcmnt_=.true.
         pposnam_=0
         if(string.ne.tab)pchar=lprefx+1
         return
         end
!
!
!
!
!
!
!
! >>>>>> Put a text sequence into the CIF.
!
         function ptext_(name,string)      
!
         logical    ptext_
         integer    lastnb
         include   'ciftbx.sys'
         logical    flag,tflag
         integer    ll
         character  name*(*),temp*(NUMCHAR),string*(*),store*(NUMCHAR)
         character  temp2*(MAXBUF)
         data store/'                                '/
!
         ptext_=.true.
         flag  =.true.
         tflag =.true.
         ll=lastnb(string)
         if(ploopf.eq.'yes') ploopc=0
         ploopf='no '
         temp=name
         if(ptextf.eq.'no ')    goto 100
         if(temp.eq.store)      goto 150
         call eotext
!
100      if(name(1:1).ne.' ')   goto 110
         if(ptextf.eq.'yes')    goto 150
         goto 130
!
110      if(ploopn.ne.0)        call eoloop
         if(vcheck.eq.'no ')    goto 120
         call dcheck(name,'char',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0) &
             temp=dictag(aroot(xdchk))
         endif
         ptext_=flag
120      pchar=-1
         if(pposnam_.ne.0) pchar=pposnam_+lprefx
         call putstr(temp)
         if(.not.flag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not in dictionary')
         endif
         if(.not.tflag) then
           if(.not.tabl_) pchar=lprefx+57
           call putstr('#< not correct type')
         endif
130      ptextf='yes'
         store=temp
         if(string(1:1).eq.' '.and.ll.gt.1) then
           pchar=-1
           temp2=';'//string(2:ll)
           call putstr(temp2)
           pchar=-1
           return
         endif
         pchar=-1
         call putstr(';')
         pchar=-1
         if(string.eq.' ') return
150      pchar=-1
         call putstr(string(1:max(1,ll)))
         pchar=-1
         pposnam_=0
         pposval_=0
         pposdec_=0
         pposend_=0
         return
         end
!
!
!
!
!
!
! >>>>>> Put a loop_ data name into the CIF.
!
         function ploop_(name)      
!
         logical    ploop_
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR)
!
         ploop_=.true.
         flag  =.true.
         if(ptextf.eq.'yes')    call eotext
         if(ploopf.eq.'no ')    call eoloop
         temp=' '
         if(name(1:1).eq.' ')   goto 100
!
         if(tabl_.and.pposnam_.eq.0) then
           temp='    '//name
         else
           temp=name
         endif
         if(vcheck.eq.'no ')    goto 100
         call dcheck(name,'    ',flag,tflag)
         if (aliaso_.and.xdchk.ne.0) then
           if (aroot(xdchk).ne.0) then
             if(tabl_.and.pposnam_.eq.0) then
               temp='    '//dictag(aroot(xdchk))
             else
               temp=dictag(aroot(xdchk))
             endif
           endif
         endif
         ploop_=flag
100      if(ploopn.ne.0)        goto 120
         ploopf='yes'
         pchar=-1
         if(pposval_.ne.0) then
           pchar=lprefx+1
           call putstr(' ')
           pchar=pposval_+lprefx
         else
           if(pposnam_.ne.0) then
             pchar=lprefx+1
             call putstr(' ')
             pchar=pposnam_+lprefx+1
           endif
         endif
         call putstr('loop_')
         pchar=-1
         if(name(1:1).eq.' ') then
           ploopn=-1
           return
         endif
120      if(pposnam_.ne.0) pchar=pposnam_+lprefx
         call putstr(temp)
         if(flag)               goto 130
         if(.not.tabl_) pchar=lprefx+57
         call putstr('#< not in dictionary')
130      pchar=lprefx+1
         ploopn=max(ploopn,0)+1
!
150      return
         end
!
!
!
!
!
! >>>>>> Create or clear a prefix string
!        Any change in the length of the prefix string flushes
!        pending text, if any,  loops and partial output lines
!
         function prefx_(strg,lstrg)
!
         logical    prefx_
         include   'ciftbx.sys'
         character  strg*(*)
         integer    lstrg,mxline
!
         mxline=MAXBUF
         if(line_.gt.0) mxline=min(line_,MAXBUF)
         if(lstrg.ne.lprefx.and.pcharl.gt.0) then
           pchar=-1
           call putstr(' ')
         endif
         if (lstrg.le.0) then
           prefx=' '
           if(pchar.ge.lprefx+1)pchar=pchar-lprefx
           lprefx=0
         else
           if(lstrg.gt.mxline) then
             call warn(' Prefix string truncated')
           endif
           prefx=strg
           if(pchar.ge.lprefx+1)pchar=pchar-lprefx+lstrg
           obuf(1:min(mxline,lstrg))=prefx
           lprefx=lstrg
           if(mxline-lprefx.lt.NUMCHAR) then
             call warn(' Output prefix may force line overflow')
           endif
         endif
         prefx_=.true.
         return
         end
!
!
!
!
!
!
! >>>>>> Close the CIF
!
         subroutine close_
!
         include   'ciftbx.sys'
!
         if(ptextf.eq.'yes') call eotext
         if(ploopn.ne.0)     call eoloop
         if(pcharl.ge.lprefx+1) then
           pchar=-1
           call putstr(' ')
         endif
         if (file_(1:1) .ne. ' ') then
           close(outdev)
           precn_=0
         endif
         return
         end
!
!
!
!
!
! >>>>>> Put the string into the output CIF buffer 
!
         subroutine putstr(string)     
!
         integer    lastnb
         include   'ciftbx.sys'
         SAVE
         character  string*(*),temp*(MAXBUF),bfill*(MAXBUF)
         character  temp2*(MAXBUF)
         integer    i,ii,mxline,ioffst,ifree,icpos,itpos
         integer    ixpos,ixtpos,it,im,kbin,kpass
         logical    pflush,waslop
         data       waslop /.false./

!
         bfill = ' '
         mxline=MAXBUF
         if(line_.gt.0) mxline=min(line_,MAXBUF)
         temp=string
         temp2=temp
         pflush=.false.
         if(pchar.lt.0) pflush=.true.
!
         do 100 i=MAXBUF,1,-1
         if(temp(i:i).eq.' ')              goto 100
         if(ptabx_.and.temp(i:i).eq.tab)    goto 100
         goto 110
100      continue
         i=0
         it=i
!
!....... Organise the output of loop_ items
!
110      if(i.eq.0)             goto 130
         if(i.eq.1.and.string.eq.tab) goto 130
         if(i.eq.1.and.string.eq.char(0)) then
           pcharl=MAXBUF
           goto 200
         endif           
         if(temp(1:1).eq.'#')   goto 130
         if(ploopf.eq.'yes')    goto 130
         if(ptextf.eq.'yes')    goto 130
         if(ploopn.le.0)        goto 130
         ploopc=ploopc+1
         if(align_.or.tabl_) then
           if(ploopc.gt.ploopn) then
             if(pcharl.gt.lprefx) pflush=.true.
             ploopc=1
             if(pchar.gt.0) pchar=lprefx+1
           endif
           if(pchar.lt.0)    goto 130
           if(tabl_) then
           kbin=(mxline-lprefx)/8
           if(ploopn.lt.kbin) then
             if(kbin/(ploopn+1).gt.1) then
             pchar=9+lprefx+ &
               (ploopc-1)*8*(kbin/(ploopn+1))
             else
             pchar=1+lprefx+ &
               (ploopc-1)*8*(kbin/ploopn)
             endif
           else
             if(ploopc.le.kbin) then
               pchar=1+lprefx+(ploopc-1)*8
             else
               kpass=(ploopc-kbin-1)/(kbin-1)+1
               pchar=2*kpass+1+lprefx+ &
                 mod(ploopc-kbin-1,kbin-1)*8
             endif
           endif
           else
             if(ptabx_) then
             icpos=1
             itpos=1
120          ixpos = 0
             if (icpos.le.i) ixpos=index(temp(icpos:i),tab)
             ixtpos=(pchar+itpos-1+ixpos)
             ixtpos=((ixtpos+7)/8)*8
             if(ixpos.gt.0) then
               if(ixpos.gt.1) then
                 temp2(itpos:ixtpos-pchar+1)=temp(icpos:ixpos-1)
               else
                 temp2(itpos:ixtpos-pchar+1)=' '
               endif
               icpos=ixpos+1
               itpos=ixtpos+2-pchar
               if(icpos.le.i) goto 120
               it=itpos-1
             else
               if(icpos.le.i) then
                 temp2(itpos:itpos+i-icpos)=temp(icpos:i)
                 it=itpos+i-icpos
               endif
             endif
             endif
             if((pchar+i).gt.mxline+1.or. &
                (ptabx_.and.pchar+it.gt.mxline+1)) then
               if(pcharl.gt.lprefx)pflush=.true.
               pchar=lprefx+1
             endif
           endif
         else
           if(ploopc.le.ploopn)   goto 130
           ploopc=1
         endif
!
!....... Is the buffer full and needs flushing?
!
130      if(i.eq.1.and.string.eq.tab) then
           if(pcharl.gt.lprefx) then
             if(obuf(pcharl:pcharl).eq.' ') pcharl=pcharl-1
           endif
         endif
         if(pchar.le.pcharl.and.pcharl.gt.lprefx) pflush=.true.
         pchar=max(lprefx+1,pchar)
         if((ploopf.eq.'yes'.or.ploopn.le.0).and.tabl_) &
           pchar=((pchar-lprefx+6)/8)*8+1+lprefx
         if(ptabx_) then
         icpos=1
         itpos=1
135      ixpos=0
         if(icpos.le.i) ixpos=index(temp(icpos:i),tab)
         ixtpos=(pchar+itpos-1+ixpos)
         ixtpos=((ixtpos+7)/8)*8
         if(ixpos.gt.0) then
           if(ixpos.gt.1) then
             temp2(itpos:ixtpos-pchar+1)=temp(icpos:ixpos-1)
           else
             temp2(itpos:ixtpos-pchar+1)=' '
           endif
           icpos=ixpos+1
           itpos=ixtpos+2-pchar
           if(icpos.le.i) goto 135
           it=itpos-1
         else
           if(icpos.le.i) then
             temp2(itpos:itpos+i-icpos)=temp(icpos:i)
             it=itpos+i-icpos
           endif
         endif
         endif
         if((pchar+i).gt.mxline+1.or. &
           (ptabx_.and.pchar+it.gt.mxline+1)) then
            pflush=.true.
            pchar=mxline+1-i
            pchar=max(lprefx+1,pchar)
         endif
         if(.not.pflush)  goto 150
140      if(pcharl.gt.lprefx) then
           if(waslop.or.(.not.tabl_)) goto 145
           ioffst=0
           pcharl=max(lastnb(obuf(1:pcharl)),lprefx+1)
           ifree=mxline-pcharl
           if(ifree.gt.0) then
           im=numtab+2
           if(numtab.gt.0.and.numtab.le.MAXTAB) then
             if(obuf(itabp(numtab):itabp(numtab)).eq.'#') &
               im=im-1
           endif
           if(ifree.ge.16.and.im.lt.4.and. &
             (obuf(1+lprefx:1+lprefx).ne.'#' &
              .and.obuf(1+lprefx:1+lprefx).ne.';' &
              .and.obuf(1+lprefx:1+lprefx).ne.'_' &
              .and.obuf(1+lprefx:1+lprefx).ne.' ' &
              .and.obuf(1+lprefx:5+lprefx).ne.'data_' &
              .and.obuf(1+lprefx:5+lprefx).ne.'save_' &
              .and.obuf(1+lprefx:5).ne.'loop_')) then
             temp(1+lprefx:pcharl)=obuf(1+lprefx:pcharl)
             obuf(1+lprefx:pcharl+8)= &
               bfill(1:8)//temp(1+lprefx:pcharl)
             ioffst = 8
             ifree=ifree-8
             pcharl=pcharl+8
           endif
           do ii=1,min(MAXTAB,numtab)
             icpos=itabp(ii)+ioffst
             if(icpos.gt.pcharl)   goto 145
             if(im.lt.4) then
             itpos=(max(icpos-lprefx, &
               ii*(mxline-lprefx)/im)+6)/8
             itpos=itpos*8+1+lprefx
             else
             itpos=(max(icpos-lprefx, &
               ii*(mxline-lprefx)/im)+4)/6
             itpos=itpos*6+1+lprefx
             endif
             if((obuf(icpos:icpos).eq.''''.or. &
                obuf(icpos:icpos).eq.'"').and. &
                itpos.gt.icpos) itpos=itpos-1
             if(itpos-icpos.gt.ifree) itpos=icpos+ifree
             if(itpos.gt.icpos) then
               temp(1:pcharl-icpos+1)= &
                 obuf(icpos:pcharl)
               if(i.lt.numtab) then
                 ixpos=itabp(ii+1)+ioffst
                 if(ixpos.gt.icpos+itpos-icpos+1) then
                   if(obuf(ixpos-(itpos-icpos+1):ixpos-1).eq. &
                     bfill(1:itpos-icpos+1)) then
                     temp(ixpos-itpos+1:pcharl-itpos+1)= &
                     obuf(ixpos:pcharl)
                     pcharl=pcharl-(itpos-icpos)
                   endif
                 endif
               endif
               obuf(icpos:pcharl+itpos-icpos)= &
                 bfill(1:itpos-icpos)//temp(1:pcharl-icpos+1)
               ifree=ifree-(itpos-icpos)
               ioffst=ioffst+itpos-icpos
               pcharl=pcharl+itpos-icpos
             endif
             if(ifree.le.0)      goto 145
           enddo
           endif
145        pcharl=max(1,lastnb(obuf))
           write(outdev,'(a)') obuf(1:pcharl)
         else
           if(precn_.gt.0) then
           if(lprefx.gt.0) then
           write(outdev,'(a)') obuf(1:lprefx)
           else
           write(outdev,'(a)')
           endif
           else
           precn_=precn_-1
           endif
         endif
         waslop=.false.
         precn_=precn_+1
         do ii = 1,MAXTAB
           itabp(ii)=0
         enddo
         numtab=0
         if(lprefx.gt.0) then
           obuf=prefx(1:lprefx)
         else
           obuf=' '
         endif
!
!....... Load the next item into the buffer
!
150      pcharl=pchar+i
         if(ptabx_) pcharl=pchar+it
         waslop= ploopf.eq.'no '.and.ploopn.gt.0.and.align_
         if(i.eq.0) then
           if(pcharl.eq.lprefx+1.and. &
             obuf(lprefx+1:lprefx+1).eq.' ') pcharl=pcharl-1
             pchar=pcharl+1
           goto 200
         endif
         if(ptabx_) then
           obuf(pchar:pcharl)=temp2(1:it)
         else
           if(string.eq.tab) pcharl=pcharl-1
           obuf(pchar:pcharl)=string(1:i)
         endif
         if(pchar.gt.1+lprefx) then
           numtab=numtab+1
           if(numtab.le.MAXTAB) itabp(numtab)=pchar
         endif
         pchar=pcharl+1
         if(pchar.gt.mxline+2) then
           call warn(' Output CIF line longer than line_')
         endif
!
200      return
         end
!
!
!
!
!
! >>>>>> Convert the number and esd to string nnnn(m), limited
!        by relative precision prec
!
         subroutine putnum(numb,sdev,prec)  
!
         include   'ciftbx.sys'
         character  string*30,temp*30,c*1,sfmt*8
         double precision numb,sdev,prec,xxnumb,xsdev,slog
         integer    i,iexp,ifp,ii,jj,j,jlnz,jn,kexp,m,ixsdev,islog
         integer    kdecp,ibexp,lexp
!
         kdecp=0
         if (sdev.gt.abs(numb)*prec) then
           if (iabs(esdlim_).ne.esdcac) then
!
!            determine the number of digits set by esdlim_
!
             if (iabs(esdlim_).lt.9 .or.iabs(esdlim_).gt.99999) then
               call warn(' Invalid value of esdlim_ reset to 19')
               esdlim_ = 19
             endif
!
!            determine the number of esd digits
!
             esddigx = 1.+alog10(float(iabs(esdlim_)))
             esdcac = iabs(esdlim_)
           endif
!
!          if esdlim_ < 0, validate pesddig_
!
           if (esdlim_.lt. 0 )then
             if (pesddig_.lt.0 .or. pesddig_.gt.5) then
               call warn(' Invalid value of pesddig_ reset to 0')
               pesddig_ = 0
             endif
           endif
!
!          determine kexp, the power of 10 necessary
!          to present sdev as an integer in the range
!          (esdlim_/10,esdlim_] or [1,-esdlim_] if esdlim_ < 0
!
           slog = dlog10(sdev)
           islog = slog+1000.
           islog = islog-1000
           kexp = -islog+esddigx
!
!          Adjust exponent kexp, so that sdev*10**kexp
!          is in the interval (esdlim_/10,esdlim_] or [1,-esdlim_]
!
 20        if (kexp.lt.minexp) then
             call warn(' Underflow of esd')
             ixsdev = 0
             go to 30
           endif
           if (kexp.gt.-minexp) then
             call warn(' Overflow of esd')
             ixsdev = 99999
             go to 30
           endif
           xsdev = sdev*10.D0**kexp
           ixsdev = xsdev+.5
           if (ixsdev.gt.iabs(esdlim_)) then
             kexp = kexp -1
             go to 20
           endif
           if (ixsdev.lt.(iabs(esdlim_)+5)/10) then
             kexp = kexp+1
             go to 20
           endif
!
!          lexp holds the number of trailing zeros which may be
!          sacrificed in the esd if the number itself has 
!          trailing zeros in the fraction which is permitted if
!          esdlim_ is negative
!
!          If esdlim_ is negative and pesddig_ is .gt.0,
!          pesddig_ will be used to force the number of digits
!          in which case lexp has the number of digits that
!          must be sacrificed (lexp > 0) or zeros to add (lexp < 0)
!           
           lexp=0
           if(esdlim_.lt.0) then
             if(pesddig_.gt.0) then
25             continue
               if(ixsdev*10**(-lexp).ge.10**(pesddig_))then
                 if(lexp.gt.0) &
                   ixsdev=ixsdev-5*10**(lexp-1)
                 ixsdev=ixsdev+5*10**lexp
                 lexp=lexp+1
                 goto 25
               endif
               if(ixsdev.lt.10**(pesddig_-1+lexp) &
                 .and.lexp.gt.0) then
                 if(ixsdev*10**(-lexp).le.iabs(esdlim_))then
                   lexp =lexp-1
                   if(lexp.ge.0) then
                     ixsdev=ixsdev-5*10**lexp
                   endif
                   if(lexp.gt.0) then
                     ixsdev=ixsdev+5*10**(lexp-1)
                   endif
                   goto 25
                 endif
               endif
               kexp=kexp-lexp
               ixsdev = ixsdev/(10**lexp)
               lexp=0
             else
             do ii = 1,4
               if(mod(ixsdev,10**ii).ne.0) go to 30
               lexp = ii
             enddo
             endif
           endif
!
!          We need to present the number to the same scaling
!          at first, but will adjust to avoid Ennn notation
!          if possible
!
 30        xxnumb = dabs(numb)*10.d0**kexp+.5
           if(xxnumb*prec .gt.1.D0) then
             call warn(' ESD less than precision of machine')
             ixsdev=0
           endif
           if(numb.lt.0.d0) xxnumb = -xxnumb
           write(string,ndpfmt)xxnumb
           if(xxnumb.lt.1.d0 .and. xxnumb.ge.0.d0)   &
              string='                         0.0E0'
           if(xxnumb.gt.-1.d0 .and. xxnumb.lt.0.d0)  &
              string='                        -0.0E0'
!
!          Extract the power of 10
!
           iexp = 0
           ibexp = 0
           do ii = 0,4
             i = 30-ii
             c = string(i:i)
             m = index('0123456789',c)
             if (m.gt.0) then
               iexp = iexp+(m-1)*10**(ii-ibexp)
             else
               if (c.eq.' ') then
                 ibexp = ibexp+1
               else
               if (c.eq.'-') iexp=-iexp
               goto 40
               endif
             endif
           enddo
           call err(' Internal error in putnum')
!
!          Scan the rest of the string shifting the
!          decimal point to get an integer
!
40         ifp = 0
           j=1
           do ii = 1,i-1
           c = string(ii:ii)
           if (c.ne.' ')then
             m=index('0123456789+-',c)
             if(m.ne.0) then
               temp(j:j)=c
               if(j.gt.1.or.c.ne.'0')j=j+1
               if(j.eq.3.and.temp(1:2).eq.'-0')j=j-1
               if(ifp.ne.0)then
                 iexp=iexp-1
                 if(iexp.le.0) goto 50
               endif
             else
               if(c.eq.'.') then
                 ifp=1
                 if(iexp.le.0) goto 50
               endif
             endif
           endif
           enddo
!
!          The string from 1 to j-1 has an integer
!          If iexp < 0, we present a 0.  If iexp > 0
!          we pad with zeros
!
50         if(j.eq.2 .and. temp(1:1).eq.'-') then
              temp(1:2)='-0'
              j=3
              iexp=0
           endif
           if(j.eq.1 .or. iexp.lt.0) then
             temp(1:1)='0'
             j=2
             iexp = 0
             if(xxnumb.lt.0.d0) then
               temp(1:2)='-0'
               j=3
             endif
           endif
           if (iexp.gt.0) then
             do ii = 1,iexp
             temp(j:j)='0'
             j=j+1
             enddo
             iexp=0
           endif
           string=temp(1:j-1)
!
!          We have the number for which the presentation
!          would be nnnnnE-kexp.  If kexp is gt 0, we can
!          decrease it and introduce a decimal point
!
           jj=0
           if(index('0123456789',temp(1:1)).eq.0) jj=1
           if(kexp.gt.0.and.kexp.lt.j-jj+8) then
             if(kexp.lt.j-1) then
               if(plzero_ .and. &
                 j-1-kexp.eq.1.and.temp(1:1).eq.'-') then
                 string=temp(1:j-1-kexp)//'0.'// &
                   temp(j-kexp:j-1)
                 j=j+2
               else
                 string=temp(1:j-1-kexp)//'.'// &
                 temp(j-kexp:j-1)
                 j=j+1
               endif
               kexp = 0
             else
               if(jj.ne.0)string(1:1)=temp(1:1)
               if(plzero_) then
                 string(1+jj:2+jj)='0.'
                 do ii=1,kexp-(j-1-jj)
                   string(2+jj+ii:2+jj+ii)='0'
                 enddo
                 string(3+jj+(kexp-(j-1-jj)):30)= &
                   temp(1+jj:j-1)
                 j=j+2+kexp-(j-1-jj)
               else
                 string(1+jj:1+jj)='.'
                 do ii=1,kexp-(j-1-jj)
                   string(1+jj+ii:1+jj+ii)='0'
                 enddo
                 string(2+jj+(kexp-(j-1-jj)):30)= &
                   temp(1+jj:j-1)
                 j=j+1+kexp-(j-1-jj)
               endif
               kexp=0
             endif
           endif
           kdecp=index(string(1:j-1),'.')
           if(kdecp.gt.0.and.kdecp.lt.j-1.and.lexp.gt.0) then
             jj=0
             do ii = 1,min(lexp,j-1-kdecp)
               c = string(j-ii:j-ii)
               if(c.ne.'0') goto 60
               jj=jj+1
             enddo
60           j=j-jj
             ixsdev=ixsdev/10**jj
             if(.not.pdecp_.and.string(j-1:j-1).eq.'.') then
               j=j-1
               kdecp=0
             endif
           endif
           if(kdecp.eq.0) then
             kdecp=j
             if(pdecp_) then
               if(plzero_.and. &
                 (j.eq.1 .or. (j.eq.2.and.string(1:1).eq.'-'))) then
                 string(j:j)='0'
                 j=j+1
               endif
               string(j:j)='.'
               j=j+1
             endif
           endif
           if(kexp.ne.0) then
             write(temp(1:5),'(i5)') -kexp
             string(j:j)='E'
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
           endif
!
!          if there is a standard deviation
!          append it in parentheses
!
           if(ixsdev.ne.0) then
             write(temp(1:5),'(i5)') ixsdev
             string(j:j)='('
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
             string(j:j)=')'
             j=j+1
           endif
         else
!
!          There is no standard deviation, just write numb
!          But limit to the digits implied by prec
!
           slog = dlog10(min(.1D0,max(prec,dpprc)))
           islog = slog+1000.5
           islog = islog-1000
           kexp = -islog
           write(sfmt,'(5h(D30.,i2,1h))') kexp
           write(temp,sfmt)numb
!
!          Now have the number in the form 
!          [sign][0].nnnnnnnnDeee
!          which, while sufficient, is not neat
!          we reformat for the case 0<=eee<=kexp
!
!
!          Extract the power of 10
!
           iexp = 0
           ibexp = 0
           do ii = 0,4
             i = 30-ii
             c = temp(i:i)
             m = index('0123456789',c)
             if (m.gt.0) then
               iexp = iexp+(m-1)*10**(ii-ibexp)
             else
               if (c.eq.' ') then
                 ibexp = ibexp+1
               else
               if (c.eq.'-') iexp=-iexp
               goto 140
               endif
             endif
           enddo
           call err(' Internal error in putnum')
!
!          Scan the rest of the string shifting the
!          decimal point to get a number with exponent 0,
!          if possible
!
140        ifp = 0
           j=1
           do ii = 1,i-1
           jn=ii
           c = temp(ii:ii)
           if (c.ne.' ')then
             m=index('0123456789+-',c)
             if(m.ne.0) then
               string(j:j)=c
               if(j.gt.1.or.c.ne.'0')j=j+1
               if(j.eq.3.and.string(1:2).eq.'-0')j=j-1
               if(ifp.ne.0)then
                 iexp=iexp-1
                 if(iexp.le.0) goto 150
               endif
             else
               if(c.eq.'.') then
                 ifp = -1
                 if(iexp.le.0) goto 150
               endif
             endif
           endif
           enddo
150        if(plzero_ .and. &
             (j.eq.1 .or.(j.eq.2.and.string(1:1).eq.'-'))) then
             string(j:j)='0'
             j=j+1
           endif
           string(j:j)='.'
           ifp = j
           j = j+1
           jlnz = j-1
155        do ii = jn+1,i-1
             c = temp(ii:ii)
             if (c.ne.' ')then
               m=index('0123456789',c)
               if(m.ne.0) then
                 string(j:j)=c
                 j=j+1
                 if(m.ne.1)jlnz=j
                 if(m.eq.1.and.ifp.ge.1.and. &
                   pposdec_.ne.0.and.pposend_.ne.0) then
                   if(j-1-ifp-min(iexp,0).le.pposend_-pposdec_) &
                     jlnz=j
                 endif
               else
                 goto 160
               endif
             endif
           enddo
160        j=jlnz
           if(j.eq.1) then
            string(1:1)='0'
            j=2
           endif
           if(iexp.lt.0.and.iexp.gt.-7.and.ifp.lt.j-1.and. &
             ifp.ne.0.and.j-ifp-iexp.le.kexp) then
             temp(1:ifp)=string(1:ifp)
             do ii = 1,-iexp
               temp(ifp+ii:ifp+ii) = '0'
             enddo
             temp(ifp-iexp+1:j-iexp-1) = string(ifp+1:j-1)
             j = j-iexp
             iexp=0
             string(1:j-1) = temp(1:j-1)
           endif
           kdecp=index(string(1:j-1),'.')
           if(kdecp.eq.0) then
             kdecp=j
             if(pdecp_) then
               string(kdecp:kdecp)='.'
               j=j+1
             endif
           endif
           if(iexp.ne.0) then
             write(temp(1:5),'(i5)')iexp
             string(j:j)='E'
             j=j+1
             do ii=1,5
               c=temp(ii:ii)
               if(c.ne.' ') then
                 string(j:j)=c
                 j=j+1
               endif
             enddo
           endif
         endif
!
         if(j.lt.1) then
           string(1:1)='0'
           j=2
         endif
         if(kdecp.lt.1)kdecp=j
         if(pposdec_.ne.0) then
           pchar=lprefx+pposdec_-kdecp+1
         else
           if(pposval_.ne.0)pchar=lprefx+pposval_
         endif
         call putstr(string(1:j-1))
         return
         end
!
!
!
!
!
! >>>>>> Check dictionary for data name validation    
!
         subroutine dcheck(name,type,flag,tflag)
!
         include   'ciftbx.sys'
         logical    flag,tflag
         character  name*(*),temp*(NUMCHAR), &
                    locase*(MAXBUF),type*4
!
         flag=.true.
         tflag=.true.
         temp=locase(name)
         call hash_find(temp, &
           dicnam,dicchain,NUMDICT,ndict,dichash,NUMHASH,xdchk)
         if(xdchk.eq.0) goto 150
         if(tcheck.eq.'no ')          goto 200
         if(type.eq.dictyp(xdchk))    goto 200
         if(type.eq.'    ')           goto 200
         if(dictyp(xdchk).eq.'text' .and. type.eq.'char') goto 200
         tflag=.false.
         goto 200
150      flag=.false.
200      continue
         return
         end
!
!
!
!
!
! >>>>>> End of text string
!
         subroutine eotext
!
         include   'ciftbx.sys'
!
         if(ptextf.ne.'yes') then
           call warn(' Out-of-sequence call to end text block')
           return
         endif
         ptextf='no '
         pchar=-1
         call putstr(';')
         call putstr(char(0))
         return
         end
!
!
!
!
!
! >>>>>> End of loop detected; check integrity and tidy up pointers
!
         subroutine eoloop
!
         include   'ciftbx.sys'
         integer   i
!
         if(ploopn.eq.0)          goto 200
         if(ploopn.eq.-1) then
           call putstr('_DUMMY')
           ploopn=1
           ploopc=0
           call warn( &
             ' Missing: missing loop_ name set as _DUMMY')
         endif
         if(ploopn.eq.ploopc)     goto 200
         do 150 i=ploopc+1,ploopn
150      call putstr('DUMMY')
         call warn(     &
               ' Missing: missing loop_ items set as DUMMY')
!
200      ploopc=0
         ploopn=0
         return
         end
!
!
!
!
!
!
! >>>>>> Set common default values
!
         block data
!
         include   'ciftbx.sys'
         data cifdev     /1/
         data outdev     /2/
         data dirdev     /3/
         data errdev     /6/
         data recbeg_    /1/
         data recend_    /0/
         data loopct     /0/
         data nhash      /0/
         data ndict      /0/
         data nname      /0/
         data nbloc      /0/
         data ploopn     /0/
         data ploopc     /0/
         data ploopf     /'no '/
         data ptextf     /'no '/
         data pfilef     /'no '/
         data testfl     /'no '/
         data vcheck     /'no '/
         data tcheck     /'no '/
         data catchk     /'yes'/
         data align_     /.true./
         data append_    /.false./
         data tabl_      /.true./
         data tabx_      /.true./
         data ptabx_     /.true./ 
         data text_      /.false./
         data loop_      /.false./
         data ndcname    /0/
         data ncname     /0/
         data save_      /.false./
         data saveo_     /.false./
         data psaveo     /.false./
         data glob_      /.false./
         data globo_     /.false./
         data alias_     /.true./
         data aliaso_    /.false./
         data nblank_    /.false./
         data nblanko_   /.false./
         data decp_      /.false./
         data pdecp_     /.false./
         data lzero_     /.false./
         data plzero_    /.false./
         data dchash     /NUMHASH*0/
         data dichash    /NUMHASH*0/
         data dhash      /NUMHASH*0/
         data dcchain    /NUMDICT*0/
         data aroot      /NUMDICT*0/
         data keychain   /NUMDICT*0/
         data ccatkey    /NUMDICT*0/
         data cindex     /NUMBLOCK*0/
         data line_      /200/
         data lastch     /0/
         data dictype_   /' '/
         data dicname_   /' '/
         data dicver_    /' '/
         data diccat_    /' '/
         data tagname_   /' '/
         data prefx      /' '/
         data tbxver_    /'CIFtbx version 2.6.2 16 Jun 1998'/
         data lprefx     /0/
         data esdlim_    /19/
         data esddig_    /0/
         data pesddig_   /0/
         data esdcac     /19/
         data esddigx    /2/
         data esdfmt     /'(e12.2)'/
         data edpfmt     /'(d12.2)'/
         data ndpfmt     /'(d30.14)'/
         data decprc     /1.e-6/
         data dpprc      /1.d-14/
         data decmin     /1.e-37/
         data dpmin      /1.d-307/
         data minexp     /-307/
         data itabp      /MAXTAB*0/
         data jrect      /-1/
         data numtab     /0/
         data recn_      /0/
         data precn_     /0/
         data posnam_    /0/
         data posval_    /0/
         data posdec_    /0/
         data posend_    /0/
         data pposnam_   /0/
         data pposval_   /0/
         data pposdec_   /0/
         data pposend_   /0/
         data quote_     /' '/
         data pquote_    /' '/
         data ibkmrk     /MAXBOOK*-1,MAXBOOK*-1, &
                          MAXBOOK*-1,MAXBOOK*-1/

         end
!
!
!       change the following include to include 'clearfp_sun.f'
!       for use on a SUN
!
        include 'clearfp.f'

