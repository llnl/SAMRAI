/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Concrete subclass of tbox::Serializable.
 *
 ************************************************************************/

#include "MainRestartData.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */

MainRestartData::MainRestartData(
   const std::string& object_name,
   std::shared_ptr<tbox::Database> input_db):
   d_object_name(object_name)
{
   TBOX_ASSERT(input_db);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   /* if not starting from restart file, set loop_time and iteration_number */
   if (!is_from_restart) {
      d_loop_time = d_start_time;
      d_iteration_number = 0;
   }
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

MainRestartData::~MainRestartData()
{
}

/*
 *************************************************************************
 *
 * Accessor methods for data
 *
 *************************************************************************
 */
int MainRestartData::getMaxTimesteps()
{
   return d_max_timesteps;
}

double MainRestartData::getStartTime()
{
   return d_start_time;
}

double MainRestartData::getEndTime()
{
   return d_end_time;
}

int MainRestartData::getRegridStep()
{
   return d_regrid_step;
}

int MainRestartData::getTagBuffer()
{
   return d_tag_buffer;
}

double MainRestartData::getLoopTime()
{
   return d_loop_time;
}

void MainRestartData::setLoopTime(
   const double loop_time)
{
   d_loop_time = loop_time;
}

int MainRestartData::getIterationNumber()
{
   return d_iteration_number;
}

void MainRestartData::setIterationNumber(
   const int iter_num)
{
   d_iteration_number = iter_num;
}

/*
 *************************************************************************
 *
 * tbox::Database input/output methods
 *
 *************************************************************************
 */
void MainRestartData::putToRestart(
   const std::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("d_max_timesteps", d_max_timesteps);
   restart_db->putDouble("d_start_time", d_start_time);
   restart_db->putDouble("d_end_time", d_end_time);
   restart_db->putInteger("d_regrid_step", d_regrid_step);
   restart_db->putInteger("d_tag_buffer", d_tag_buffer);
   restart_db->putDouble("d_loop_time", d_loop_time);
   restart_db->putInteger("d_iteration_number", d_iteration_number);
}

void MainRestartData::getFromInput(
   std::shared_ptr<tbox::Database> input_db,
   bool is_from_restart)
{
   TBOX_ASSERT(input_db);

   if (input_db->keyExists("max_timesteps")) {
      d_max_timesteps = input_db->getInteger("max_timesteps");
   } else {
      if (!is_from_restart) {
         TBOX_ERROR("max_timesteps not entered in input file" << std::endl);
      }
   }

   if (input_db->keyExists("start_time")) {
      d_start_time = input_db->getDouble("start_time");
   } else {
      d_start_time = 0.0;
   }

   if (input_db->keyExists("end_time")) {
      d_end_time = input_db->getDouble("end_time");
   } else {
      d_end_time = 100000.;
   }

   if (input_db->keyExists("regrid_step")) {
      d_regrid_step = input_db->getInteger("regrid_step");
   } else {
      d_regrid_step = 2;
   }

   if (input_db->keyExists("tag_buffer")) {
      d_tag_buffer = input_db->getInteger("tag_buffer");
   } else {
      d_tag_buffer = d_regrid_step;
   }
}

void MainRestartData::getFromRestart()
{
   std::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in the restart file.");
   }
   std::shared_ptr<tbox::Database> restart_db(
      root_db->getDatabase(d_object_name));

   d_max_timesteps = restart_db->getInteger("d_max_timesteps");
   d_start_time = restart_db->getDouble("d_start_time");
   d_end_time = restart_db->getDouble("d_end_time");
   d_regrid_step = restart_db->getInteger("d_regrid_step");
   d_tag_buffer = restart_db->getInteger("d_tag_buffer");
   d_loop_time = restart_db->getDouble("d_loop_time");
   d_iteration_number = restart_db->getInteger("d_iteration_number");
}
