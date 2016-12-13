#ifndef GRAPHITE_ADJUDICATOR_MAPPING_H
#define GRAPHITE_ADJUDICATOR_MAPPING_H

#include "core/mapping/IMapping.h"
#include "core/adjudicator/IAdjudicator.h"

#include "gssw/gssw.h"

namespace graphite
{
namespace adjudicator
{

	class GSSWMapping : public IMapping
	{
    public:
		typedef std::shared_ptr< GSSWMapping > SharedPtr;
        GSSWMapping(std::shared_ptr< gssw_graph_mapping > gsswMappingPtr, IAlignment::SharedPtr alignmentPtr);
		~GSSWMapping();

		int getMappingScore() override;
		std::shared_ptr< gssw_graph_mapping > getGSSWMappingStruct();
		/* MappingAlignmentInfo::SharedPtr getMappingAlignmentInfo(IAllele::SharedPtr allelePtr, IAdjudicator::SharedPtr adjudicatorPtr) override; */
		IAlignment::SharedPtr getAlignmentPtr() override;
		std::vector< IAllele::SharedPtr > getAllelePtrs() override;
		position getPosition() override { return m_position; }
		std::vector< MappingAlignmentInfo::SharedPtr > getMappingAlignmentInfoPtrs(IAdjudicator::SharedPtr adjudicatorPtr);
		void incrementAlleleCounts() override;
		void outputSAMEntry(IAlignment::SharedPtr alignmentPtr) override;
		void setMapped(bool mapped) override;
		bool getMapped() override { return m_mapped; }
		void addAlleleCountCallback(std::function< void () > functor) override;

		void printMapping() override;
		void printSimpleMapping();

    private:

		// std::vector<std::pair<int32_t, char>> removeCigarXs(std::vector<std::pair<int32_t, char>> raw_cigar);
		std::vector<std::pair<int32_t, char>> extractCigar();
		std::string stringifySimpleCigar(std::vector<std::pair<int32_t, char>> simple_cigar);
		std::shared_ptr< gssw_graph_mapping > m_gssw_mapping_ptr;
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		std::unordered_map< IAllele*, gssw_node* > m_allele_gssw_nodes_map;
		IAlignment::SharedPtr m_alignment_ptr;
		position m_position;
		std::vector< std::function< void () > > m_allele_incrementor_callback_list;
		bool m_mapped;
	};

}
}

#endif //GRAPHITE_GSSW_MAPPING_H
