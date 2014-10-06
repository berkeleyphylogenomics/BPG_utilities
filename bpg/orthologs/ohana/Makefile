CXX = g++
CXXLD = g++
CXXFLAGS = -O2
CXXCOMPILE = $(CXX) $(CXXFLAGS)
CXXLINK = $(CXXLD) $(CXXFLAGS) -o $@
LIBS = -lboost_program_options -lpq -lpqxx -lbpp-phyl -lboost_regex

find_orthologs: main.o
	@rm -f find_orthologs
	$(CXXLINK) main.o modified_preorder_tree_traversal.o put_left_right_ids_in_db.o get_taxon_of_leaf_map.o get_breadth_first_visit_order.o get_longest_distance_in_tree.o find_orthologs_in_tree.o find_inparalogs_in_tree.o put_inparalogs_in_db.o find_proximal_subtrees.o put_orthologs_in_db.o put_superorthologous_nodes_in_db.o $(LIBS)

main.o: main.cpp modified_preorder_tree_traversal.o put_left_right_ids_in_db.o get_taxon_of_leaf_map.o get_breadth_first_visit_order.o get_longest_distance_in_tree.o find_orthologs_in_tree.o find_inparalogs_in_tree.o put_inparalogs_in_db.o find_proximal_subtrees.o put_orthologs_in_db.o put_superorthologous_nodes_in_db.o

.cpp.o:
	$(CXXCOMPILE) $(INCLUDE_PATH) -c $<

clean:
	@rm -f *.o find_orthologs

