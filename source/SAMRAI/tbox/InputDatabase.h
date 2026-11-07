/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   An input database structure that stores (key,value) pairs
 *
 ************************************************************************/

#ifndef included_tbox_InputDatabase
#define included_tbox_InputDatabase

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/MemoryDatabase.h"

namespace SAMRAI {
namespace tbox {

/**
 * @brief Class InputDatabase stores (key,value) pairs in a hierarchical
 * database.
 *
 * This is just another name for the MemoryDatabase. @see MemoryDatabase
 *
 * It is normally filled with data using a Parser (@see Parser) and used to
 * pass user supplied input from input file to constructors for problem setup.
 *
 */
typedef MemoryDatabase InputDatabase;

}
}

#endif
