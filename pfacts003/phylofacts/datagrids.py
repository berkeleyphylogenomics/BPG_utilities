"""Keyword Search DataGrids Module

    The pfacts003.datagrids.datagrids module provides BPG Columns and the
Related Data Grid. This module inherits from RelatedDataGrid to provide the
FamilyDataGrid.
"""

from pfacts003.utils.datagrids.datagrids import *
from djblets.datagrid.grids import DataGrid


class FamilyDataGrid(RelatedDataGrid):
    family_id = Column(label="Family ID", field_name="id", db_field="id", sortable=False)
    n_subfamilies = Column("Number of Subfamilies", sortable=False)
    max_length = Column(
        field_name="length_longest",
        label="Max Length",
        sortable=False,
    )
    scientific_name = ForeignKeyColumn(label="Species (Scientific Name)",
        field_name="most_recent_common_taxid.scientific_name",
        sortable=False)
    common_name = ForeignKeyColumn(label="Species (Common Name)",
        field_name="most_recent_common_taxid.common_name",
        sortable=False)


    def __init__(self, request, queryset=None, title="", extra_context={},
                 optimize_sorts=True):
        DataGrid.__init__(self, request, queryset, title, extra_context)
        self.default_sort = []
        self.default_columns = [
            'family_id',
            'n_subfamilies',
            'max_length', 
        ]
        self.listview_template =\
                               "datagrid/bpg_customized/ortholog_listview.html"


