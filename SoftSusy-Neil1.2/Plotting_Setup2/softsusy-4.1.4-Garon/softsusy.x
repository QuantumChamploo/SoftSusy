#! /bin/sh

# softsusy.x - temporary wrapper script for .libs/softsusy.x
# Generated by libtool (GNU libtool) 2.4.6 Debian-2.4.6-0.1
#
# The softsusy.x program cannot be directly executed until all the libtool
# libraries that it depends on are installed.
#
# This wrapper script should never be moved out of the build directory.
# If it is, it will not operate correctly.

# Sed substitution that helps us do robust quoting.  It backslashifies
# metacharacters that are still active within double-quoted strings.
sed_quote_subst='s|\([`"$\\]\)|\\\1|g'

# Be Bourne compatible
if test -n "${ZSH_VERSION+set}" && (emulate sh) >/dev/null 2>&1; then
  emulate sh
  NULLCMD=:
  # Zsh 3.x and 4.x performs word splitting on ${1+"$@"}, which
  # is contrary to our usage.  Disable this feature.
  alias -g '${1+"$@"}'='"$@"'
  setopt NO_GLOB_SUBST
else
  case `(set -o) 2>/dev/null` in *posix*) set -o posix;; esac
fi
BIN_SH=xpg4; export BIN_SH # for Tru64
DUALCASE=1; export DUALCASE # for MKS sh

# The HP-UX ksh and POSIX shell print the target directory to stdout
# if CDPATH is set.
(unset CDPATH) >/dev/null 2>&1 && unset CDPATH

relink_command=""

# This environment variable determines our operation mode.
if test "$libtool_install_magic" = "%%%MAGIC variable%%%"; then
  # install mode needs the following variables:
  generated_by_libtool_version='2.4.6'
  notinst_deplibs=' libsoft.la libsupermodel.la libtsil.la'
else
  # When we are sourced in execute mode, $file and $ECHO are already set.
  if test "$libtool_execute_magic" != "%%%MAGIC variable%%%"; then
    file="$0"

# A function that is used when there is no print builtin or printf.
func_fallback_echo ()
{
  eval 'cat <<_LTECHO_EOF
$1
_LTECHO_EOF'
}
    ECHO="printf %s\\n"
  fi

# Very basic option parsing. These options are (a) specific to
# the libtool wrapper, (b) are identical between the wrapper
# /script/ and the wrapper /executable/ that is used only on
# windows platforms, and (c) all begin with the string --lt-
# (application programs are unlikely to have options that match
# this pattern).
#
# There are only two supported options: --lt-debug and
# --lt-dump-script. There is, deliberately, no --lt-help.
#
# The first argument to this parsing function should be the
# script's ./libtool value, followed by no.
lt_option_debug=
func_parse_lt_options ()
{
  lt_script_arg0=$0
  shift
  for lt_opt
  do
    case "$lt_opt" in
    --lt-debug) lt_option_debug=1 ;;
    --lt-dump-script)
        lt_dump_D=`$ECHO "X$lt_script_arg0" | /usr/bin/sed -e 's/^X//' -e 's%/[^/]*$%%'`
        test "X$lt_dump_D" = "X$lt_script_arg0" && lt_dump_D=.
        lt_dump_F=`$ECHO "X$lt_script_arg0" | /usr/bin/sed -e 's/^X//' -e 's%^.*/%%'`
        cat "$lt_dump_D/$lt_dump_F"
        exit 0
      ;;
    --lt-*)
        $ECHO "Unrecognized --lt- option: '$lt_opt'" 1>&2
        exit 1
      ;;
    esac
  done

  # Print the debug banner immediately:
  if test -n "$lt_option_debug"; then
    echo "softsusy.x:softsusy.x:$LINENO: libtool wrapper (GNU libtool) 2.4.6 Debian-2.4.6-0.1" 1>&2
  fi
}

# Used when --lt-debug. Prints its arguments to stdout
# (redirection is the responsibility of the caller)
func_lt_dump_args ()
{
  lt_dump_args_N=1;
  for lt_arg
  do
    $ECHO "softsusy.x:softsusy.x:$LINENO: newargv[$lt_dump_args_N]: $lt_arg"
    lt_dump_args_N=`expr $lt_dump_args_N + 1`
  done
}

# Core function for launching the target application
func_exec_program_core ()
{

      if test -n "$lt_option_debug"; then
        $ECHO "softsusy.x:softsusy.x:$LINENO: newargv[0]: $progdir/$program" 1>&2
        func_lt_dump_args ${1+"$@"} 1>&2
      fi
      exec "$progdir/$program" ${1+"$@"}

      $ECHO "$0: cannot exec $program $*" 1>&2
      exit 1
}

# A function to encapsulate launching the target application
# Strips options in the --lt-* namespace from $@ and
# launches target application with the remaining arguments.
func_exec_program ()
{
  case " $* " in
  *\ --lt-*)
    for lt_wr_arg
    do
      case $lt_wr_arg in
      --lt-*) ;;
      *) set x "$@" "$lt_wr_arg"; shift;;
      esac
      shift
    done ;;
  esac
  func_exec_program_core ${1+"$@"}
}

  # Parse options
  func_parse_lt_options "$0" ${1+"$@"}

  # Find the directory that this script lives in.
  thisdir=`$ECHO "$file" | /usr/bin/sed 's%/[^/]*$%%'`
  test "x$thisdir" = "x$file" && thisdir=.

  # Follow symbolic links until we get to the real thisdir.
  file=`ls -ld "$file" | /usr/bin/sed -n 's/.*-> //p'`
  while test -n "$file"; do
    destdir=`$ECHO "$file" | /usr/bin/sed 's%/[^/]*$%%'`

    # If there was a directory component, then change thisdir.
    if test "x$destdir" != "x$file"; then
      case "$destdir" in
      [\\/]* | [A-Za-z]:[\\/]*) thisdir="$destdir" ;;
      *) thisdir="$thisdir/$destdir" ;;
      esac
    fi

    file=`$ECHO "$file" | /usr/bin/sed 's%^.*/%%'`
    file=`ls -ld "$thisdir/$file" | /usr/bin/sed -n 's/.*-> //p'`
  done

  # Usually 'no', except on cygwin/mingw when embedded into
  # the cwrapper.
  WRAPPER_SCRIPT_BELONGS_IN_OBJDIR=no
  if test "$WRAPPER_SCRIPT_BELONGS_IN_OBJDIR" = "yes"; then
    # special case for '.'
    if test "$thisdir" = "."; then
      thisdir=`pwd`
    fi
    # remove .libs from thisdir
    case "$thisdir" in
    *[\\/].libs ) thisdir=`$ECHO "$thisdir" | /usr/bin/sed 's%[\\/][^\\/]*$%%'` ;;
    .libs )   thisdir=. ;;
    esac
  fi

  # Try to get the absolute directory name.
  absdir=`cd "$thisdir" && pwd`
  test -n "$absdir" && thisdir="$absdir"

  program='softsusy.x'
  progdir="$thisdir/.libs"


  if test -f "$progdir/$program"; then
    # Add our own library path to DYLD_LIBRARY_PATH
    DYLD_LIBRARY_PATH="/Users/glyph/Downloads/SoftSusy-master/Plotting_Setup2/softsusy-4.1.4-Garon/.libs:$DYLD_LIBRARY_PATH"

    # Some systems cannot cope with colon-terminated DYLD_LIBRARY_PATH
    # The second colon is a workaround for a bug in BeOS R4 sed
    DYLD_LIBRARY_PATH=`$ECHO "$DYLD_LIBRARY_PATH" | /usr/bin/sed 's/::*$//'`

    export DYLD_LIBRARY_PATH

    if test "$libtool_execute_magic" != "%%%MAGIC variable%%%"; then
      # Run the actual program with our arguments.
      func_exec_program ${1+"$@"}
    fi
  else
    # The program doesn't exist.
    $ECHO "$0: error: '$progdir/$program' does not exist" 1>&2
    $ECHO "This script is just a wrapper for $program." 1>&2
    $ECHO "See the libtool documentation for more information." 1>&2
    exit 1
  fi
fi
