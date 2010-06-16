#!/bin/sh

# If you're me (sourceforge: ctrapnell), use this to push the
# contents of this directory to the appropriate place on sourceforge.
# Changes will be instantaneous on http://tophat-bio.sf.net.
scp -r * ctrapnell@shell.sourceforge.net:~/tophat-bio/htdocs
