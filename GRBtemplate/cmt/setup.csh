# echo "Setting GRBtemplate v1r0p2 in /data0/glast/ScienceTools/celestialSources"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /data0/glast/extlib/rh9_gcc32/CMT/v1r16p20040701
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=GRBtemplate -version=v1r0p2 -path=/data0/glast/ScienceTools/celestialSources  $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

