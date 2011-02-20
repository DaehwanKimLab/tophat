#ifndef SEQAN_STRING_JOURNAL_BASE_H
#define SEQAN_STRING_JOURNAL_BASE_H

#define NDEBUG

#ifndef NDEBUG
#  define DEBUG_OUT(args) \
   std::cerr << args << std::endl;
#else
#  define DEBUG_OUT(args)
#endif

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

#include "string_journal_forwards.h"
#include "string_journal_operation.h"
#include "string_journal_node.h"
#include "journal.h"
#include "string_journal.h"
#include "iterator_journal.h"
#include "string_journal_interface.h"
#include "string_journal_debug.h"
#include "string_journal_utility.h"
#include "string_journal_test_foundry.h"

#endif // ndef(SEQAN_STRING_JOURNAL_BASE_H)
