# echo "Setting GRBtemplate v0r1 in /data0/glast/ScienceTools_v5/celestialSources"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /data0/glast/extlib/rh9_gcc32/CMT/v1r16p20040701
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=GRBtemplate -version=v0r1 -path=/data0/glast/ScienceTools_v5/celestialSources  $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}
