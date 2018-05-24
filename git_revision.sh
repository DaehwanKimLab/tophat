#! /bin/sh
# Display the SHA1 of the commit in which configure.ac was last modified.
# If it's not checked in yet, use the SHA1 of HEAD plus -dirty.

if [ ! -d .git ] ; then
  # if no .git directory, assume they're not using Git
  printf 'unknown'
else # show the commited git version
  printf '%s' `git show -s --pretty=format:%h`
fi

