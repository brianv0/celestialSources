# echo "Setting DarkMatter v0 in /nfs/slac/g/glast/users/glground/cohen/GRBAnalysis/celestialSources"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/slac.stanford.edu/g/glast/applications/CMT/v1r16p20040701; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh

tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=DarkMatter -version=v0 -path=/nfs/slac/g/glast/users/glground/cohen/GRBAnalysis/celestialSources  $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

