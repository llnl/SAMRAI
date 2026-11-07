/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Collects and writes out data on communication graphs.
 *
 ************************************************************************/
#include "SAMRAI/tbox/CommGraphWriter.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

#include "SAMRAI/tbox/MessageStream.h"

namespace SAMRAI {
namespace tbox {

/*
 ***********************************************************************
 ***********************************************************************
 */
CommGraphWriter::CommGraphWriter():
   d_root_rank(0),
   d_write_full_graph(true)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
CommGraphWriter::~CommGraphWriter()
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
size_t CommGraphWriter::addRecord(
   const SAMRAI_MPI& mpi,
   size_t number_of_edge_types,
   size_t number_of_node_value_types)
{
   d_records.resize(1 + d_records.size());
   Record& record = d_records.back();
   record.d_mpi = mpi;
   record.d_edges.resize(number_of_edge_types);
   record.d_node_values.resize(number_of_node_value_types);
   return d_records.size() - 1;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void CommGraphWriter::setEdgeInCurrentRecord(
   size_t edge_type_index,
   const std::string& edge_label,
   double edge_value,
   EdgeDirection edge_direction,
   int other_node)
{
   TBOX_ASSERT(edge_type_index < d_records.back().d_edges.size());

   Edge& edge = d_records.back().d_edges[edge_type_index];

   edge.d_label = edge_label;
   edge.d_value = edge_value;
   edge.d_dir = edge_direction;
   edge.d_other_node = other_node;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void CommGraphWriter::setNodeValueInCurrentRecord(
   size_t nodevalue_type_index,
   const std::string& nodevalue_label,
   double node_value)
{
   TBOX_ASSERT(nodevalue_type_index < d_records.back().d_node_values.size());

   NodeValue& nodevalue = d_records.back().d_node_values[nodevalue_type_index];

   nodevalue.d_value = node_value;
   nodevalue.d_label = nodevalue_label;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void CommGraphWriter::writeGraphToTextStream(
   size_t record_number,
   std::ostream& os) const
{
   if (d_write_full_graph) {
      writeFullGraphToTextStream(record_number, os);
   } else {
      writeGraphSummaryToTextStream(record_number, os);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void CommGraphWriter::writeGraphSummaryToTextStream(
   size_t record_number,
   std::ostream& os) const
{
   /*
    * Gather graph data on d_root_rank and write out.
    */
   TBOX_ASSERT(record_number < d_records.size());

   const Record& record = d_records[record_number];

   MessageStream ostr;
   std::vector<double> values;
   values.reserve(record.d_node_values.size() + record.d_edges.size());

   for (size_t inodev = 0; inodev < record.d_node_values.size(); ++inodev) {
      values.push_back(record.d_node_values[inodev].d_value);
   }
   for (size_t iedge = 0; iedge < record.d_edges.size(); ++iedge) {
      values.push_back(record.d_edges[iedge].d_value);
   }

   if (values.size() > 0) {
      std::vector<double> tmpvalues(values);
      record.d_mpi.Reduce(
         (void *)&tmpvalues[0],
         (void *)&values[0],
         int(values.size()),
         MPI_DOUBLE,
         MPI_MAX,
         d_root_rank);
   }

   os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
   os.precision(8);

   std::vector<NodeValue> max_nodev(record.d_node_values.size());
   std::vector<Edge> max_edge(record.d_edges.size());

   if (record.d_mpi.getRank() == d_root_rank) {

      std::vector<double>::const_iterator vi = values.begin();

      os << "\nCommGraphWriter begin record number " << record_number << '\n';
      os << "Node maximums:\n";
      for (size_t inodev = 0; inodev < record.d_node_values.size(); ++inodev) {
         os << '\t' << record.d_node_values[inodev].d_label << '\t' << *(vi++) << '\n';
      }
      os << "Edge maximums:\n";
      for (size_t iedge = 0; iedge < record.d_edges.size(); ++iedge) {
         os << '\t' << record.d_edges[iedge].d_label << '\t' << *(vi++) << '\n';
      }
      os << "CommGraphWriter end record number " << record_number << '\n';

   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void CommGraphWriter::writeFullGraphToTextStream(
   size_t record_number,
   std::ostream& os) const
{
   /*
    * Gather graph data on d_root_rank and write out.
    */
   TBOX_ASSERT(record_number < d_records.size());

   const Record& record = d_records[record_number];

   MessageStream ostr;

   for (size_t inodev = 0; inodev < record.d_node_values.size(); ++inodev) {
      const NodeValue& nodev = record.d_node_values[inodev];
      ostr << nodev.d_value;
   }
   for (size_t iedge = 0; iedge < record.d_edges.size(); ++iedge) {
      const Edge& edge = record.d_edges[iedge];
      ostr << edge.d_value << edge.d_dir << edge.d_other_node;
   }

   std::vector<char> tmpbuf(record.d_mpi.getRank() == d_root_rank ?
                            ostr.getCurrentSize() * record.d_mpi.getSize() : 0);

   if (ostr.getCurrentSize() > 0) {
      record.d_mpi.Gather(
         (void *)ostr.getBufferStart(),
         int(ostr.getCurrentSize()),
         MPI_CHAR,
         (record.d_mpi.getRank() == d_root_rank ? &tmpbuf[0] : 0),
         int(record.d_mpi.getRank() == d_root_rank ? ostr.getCurrentSize() : 0),
         MPI_CHAR,
         d_root_rank);
   }

   os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
   os.precision(8);

   std::vector<NodeValue> max_nodev(record.d_node_values.size());
   std::vector<Edge> max_edge(record.d_edges.size());

   if (record.d_mpi.getRank() == d_root_rank) {

      os << "\nCommGraphWriter begin record number " << record_number << '\n';
      os << "# proc" << '\t' << "dir" << '\t' << "remote" << '\t' << "value" << '\t' << "label\n";

      if (!tmpbuf.empty()) {
         MessageStream istr(tmpbuf.size(),
                                  MessageStream::Read,
                                  &tmpbuf[0],
                                  false);

         for (int src_rank = 0; src_rank < record.d_mpi.getSize(); ++src_rank) {

            NodeValue tmpnodev;
            for (size_t inodev = 0; inodev < record.d_node_values.size(); ++inodev) {
               istr >> tmpnodev.d_value;
               os << src_rank
                  << '\t' << tmpnodev.d_value
                  << '\t' << record.d_node_values[inodev].d_label
                  << '\n';
               if (max_nodev[inodev].d_value < tmpnodev.d_value) {
                  max_nodev[inodev] = tmpnodev;
               }
            }

            Edge tmpedge;
            for (size_t iedge = 0; iedge < record.d_edges.size(); ++iedge) {
               istr >> tmpedge.d_value >> tmpedge.d_dir >> tmpedge.d_other_node;
               os << src_rank
                  << '\t' << (tmpedge.d_dir == FROM ? "<-" : "->")
                  << '\t' << tmpedge.d_other_node
                  << '\t' << tmpedge.d_value
                  << '\t' << record.d_edges[iedge].d_label
                  << '\n';
               if (max_edge[iedge].d_value < tmpedge.d_value) {
                  max_edge[iedge] = tmpedge;
               }
            }

         }
      }

      os << "Node maximums:\n";
      for (size_t inodev = 0; inodev < record.d_node_values.size(); ++inodev) {
         os << '\t' << record.d_node_values[inodev].d_label << '\t' << max_nodev[inodev].d_value
            << '\n';
      }
      os << "Edge maximums:\n";
      for (size_t iedge = 0; iedge < record.d_edges.size(); ++iedge) {
         os << '\t' << record.d_edges[iedge].d_label << '\t' << max_edge[iedge].d_value << '\n';
      }

      os << "CommGraphWriter end record number " << record_number << '\n';

   }
}

}
}
#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
