/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Routine to time some routines in the dummy class Foo.
 *
 ************************************************************************/
#include "Foo.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include <memory>

#define LOOP_MAX (10000000)

Foo::Foo()
{
   d_depth = 0;
   d_max_depth = 0;
}

Foo::~Foo()
{
}

void Foo::timerOff()
{
   std::shared_ptr<tbox::Timer> timer(tbox::TimerManager::getManager()->
                                        getTimer("dummy::SomeClassName::shouldBeTurnedOff"));
   timer->start();

   timer->stop();
}

void Foo::timerOn()
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::timerOn()"));
   timer->start();

   timer->stop();
}

void Foo::zero(
   int depth)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::zero()"));
   if (depth > 0) {
      timer->start();
      one(depth);
      timer->stop();
   }
}

void Foo::one(
   int depth)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::one()"));
   if (depth > 1) {
      timer->start();
      two(depth);
      timer->stop();
   }
}

void Foo::two(
   int depth)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::two()"));
   if (depth > 2) {
      timer->start();
      three(depth);
      timer->stop();
   }
}

void Foo::three(
   int depth)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::three()"));
   if (depth > 3) {
      timer->start();
      four(depth);
      timer->stop();
   }
}

void Foo::four(
   int depth)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::four()"));
   if (depth > 4) {
      timer->start();
      five(depth);
      timer->stop();
   }
}

void Foo::five(
   int depth)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::five()"));
   if (depth > 5) {
      timer->start();
      six(depth);
      timer->stop();
   }
}

void Foo::six(
   int depth)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::six()"));
   if (depth > 6) {
      timer->start();
      seven(depth);
      timer->stop();
   }
}

void Foo::seven(
   int depth)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer("apps::Foo::seven()"));

   NULL_USE(timer);

   if (depth > 7) {
      TBOX_ERROR("Seven levels is maximum implemented in Foo."
         << "\n please reset exclusive_time_level to something <= 7.");
   }
}

void Foo::startAndStop(
   std::string& name)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer(name));
   timer->start();

   timer->stop();
}

void Foo::setMaxDepth(
   int max_depth)
{
   d_max_depth = max_depth;
}

void Foo::start(
   std::string& name)
{
   ++d_depth;

   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer(name));
   timer->start();

#ifndef LACKS_SSTREAM
   std::ostringstream osst;
   osst << d_depth;
   std::string out = osst.str();
   if (d_depth < d_max_depth) {
      start(out);
   }

   stop(out);
#endif
}

void Foo::stop(
   std::string& name)
{
   std::shared_ptr<tbox::Timer> timer(
      tbox::TimerManager::getManager()->getTimer(name));
   timer->stop();
}
