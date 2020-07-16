# Script to add a MOSRS branch to your local git repository


if [ "$#" -lt 2 ]; then
  echo
  echo "usage: $0 remote_svn_path local_branch_name"
  echo
  echo "   remote_svn_path:    the path on MOSRS where the branch you want to work on resides, e.g."
  echo
  echo "                        trunk"
  echo "                        branches/pkg/adrianhill/vn0.9.0_to_0.9.x_part1_pkg"
  echo "                        branches/dev/toddjones/r7304_tvfff/vn0.9.0_to_0.9.x_part1_pkg"
  echo
  echo "   local_branch_name:  the name you'd like you branch to have locally, e.g."
  echo
  echo "                        mosrs-trunk"
  echo "                        adrianhill-vn0.9.0"
  echo "                        toddjones-tvfff"
  echo
  echo "NOTE: before fetching the SVN branch your local 'master' branch"
  echo "      will be made up-to-date with the 'master' branch of the"
  echo "      Leeds MONC fork on github (https://github.com/Leeds-MONC/monc)"
  exit 1
fi

REMOTE_SVN_PATH="$1"
LOCAL_BRANCH_NAME="$2"

if [ `git branch --list ${LOCAL_BRANCH_NAME}` ];
then
  echo "The local branch '${LOCAL_BRANCH_NAME}' already exists"
  exit 1
fi

# ensure local master is up-to-date with upstream
if git ls-remote --exit-code upstream &> /dev/null; then
  echo
else
  git remote add upstream https://github.com/Leeds-MONC/monc
fi

git checkout master --quiet
git pull upstream master --quiet

if [ `git config svn-remote.${LOCAL_BRANCH_NAME}.url` ]; then
  echo
  echo "Some of the configuration for this SVN remote is already set"
  echo "(maybe you've tried running this command already, but stopped it running?)"
  echo
  echo "Please delete this configuration by issuing the following commands:"
  echo "  git config --remove-section svn-remote.${LOCAL_BRANCH_NAME}"
  exit 1
fi

git config --add svn-remote.${LOCAL_BRANCH_NAME}.url https://code.metoffice.gov.uk/svn/monc/main/branches/pkg/adrianhill/vn0.9.0_to_0.9.x_part1_pkg
git config --add svn-remote.${LOCAL_BRANCH_NAME}.fetch :refs/remotes/mosrs/${LOCAL_BRANCH_NAME}
git svn fetch ${LOCAL_BRANCH_NAME}
git checkout remotes/mosrs/${LOCAL_BRANCH_NAME} -b ${LOCAL_BRANCH_NAME}
# ensure that the history syncs up so that revisions already associated with
# git commits are identified
git rebase master
