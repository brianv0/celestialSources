# echo "Setting DarkMatter v0 in /nfs/slac/g/glast/users/glground/cohen/GRBAnalysis/celestialSources"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/slac.stanford.edu/g/glast/applications/CMT/v1r16p20040701
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=DarkMatter -version=v0 -path=/nfs/slac/g/glast/users/glground/cohen/GRBAnalysis/celestialSources  $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

