/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Collects and writes out data on communication graphs.
 *
 ************************************************************************/
#ifndef included_tbox_CommGraphWriter
#define included_tbox_CommGraphWriter

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <string>
#include <vector>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Collects data on distributed communication graphs and writes
 * out for analysis.
 *
 * A node can have multiple values, each with a label.  An node can
 * have multiple edges, each with a label.
 */
class CommGraphWriter
{

public:
   /*!
    * @brief Default constructor.
    */
   CommGraphWriter();

   /*!
    * @brief Destructor.
    */
   virtual ~CommGraphWriter();

   /*!
    * @brief Set whether to write full graph.
    *
    * Writing full graph is unscalable, but can be done at large scales
    * if you have enough computing time and memory.  Writing full graph
    * is on by default.
    */
   void setWriteFullGraph(bool write_full_graph) {
      d_write_full_graph = write_full_graph;
   }

   /*!
    * @brief Add a graph record.
    *
    * @param[in] mpi Where the graph data is distributed.
    *
    * @param[in] number_of_edge_types
    *
    * @param[in] number_of_node_value_types
    *
    * @return Index of the record.
    */
   size_t
   addRecord(
      const SAMRAI_MPI& mpi,
      size_t number_of_edge_types,
      size_t number_of_node_value_types);

   /*!
    * @brief Get the current number of records.
    *
    * @return Current number of records.
    */
   size_t getNumberOfRecords() const {
      return d_records.size();
   }

   enum EdgeDirection { FROM = 0, TO = 1 };

   /*!
    * @brief Set an edge in the current record.
    *
    * The label only matters on the root process.  Other processes do
    * nothing in this method.
    */
   void
   setEdgeInCurrentRecord(
      size_t edge_type_index,
      const std::string& edge_label,
      double edge_value,
      EdgeDirection edge_direction,
      int other_node);

   /*!
    * @brief Set a node value in the current record.
    *
    * The label only matters on the root process.  Other processes do
    * nothing in this method.
    */
   void
   setNodeValueInCurrentRecord(
      size_t nodevalue_type_index,
      const std::string& nodevalue_label,
      double node_value);

   /*!
    * @brief Gather data onto the root process and write out text file.
    */
   void
   writeGraphToTextStream(
      size_t record_number,
      std::ostream& os) const;

   struct Edge {
      Edge():d_value(0.0),
         d_dir(TO),
         d_other_node(-1) {
      }
      double d_value;
      EdgeDirection d_dir;
      int d_other_node;
      std::string d_label;
   };


private:
   // Unimplemented copy constructor.
   CommGraphWriter(
      const CommGraphWriter& other);

   // Unimplemented assignment operator.
   CommGraphWriter&
   operator = (
      const CommGraphWriter& rhs);

   struct NodeValue {
      NodeValue():d_value(0.0) {
      }
      double d_value;
      std::string d_label;
   };

   struct Record {
      Record():d_mpi(MPI_COMM_NULL) {
      }
      SAMRAI_MPI d_mpi;
      std::vector<Edge> d_edges;
      std::vector<NodeValue> d_node_values;
   };

   void
   writeGraphSummaryToTextStream(
      size_t record_number,
      std::ostream& os) const;

   void
   writeFullGraphToTextStream(
      size_t record_number,
      std::ostream& os) const;

   int d_root_rank;
   std::vector<Record> d_records;

   bool d_write_full_graph;

};

}
}

#endif  // included_tbox_CommGraphWriter
