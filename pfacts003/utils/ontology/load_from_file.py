import os
os.environ["DJANGO_SETTINGS_MODULE"] = "pfacts003.settings"

from pfacts003.utils.ontology.parsers import parse_obo_xml, parse_go_annotation_file, parse_orthodb_annotation_file
from pfacts003.utils.ontology.loader import load_ontology, load_annotations, create_annotations


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser( description='''  ''')
    parser.add_argument('--obo-xml-file', metavar='file', 
        help='OBO XML File', type=argparse.FileType('r') )
    parser.add_argument('--go-annotation-file', metavar='file', 
        help='File in GO Annotation Format', type=argparse.FileType('r') )
    parser.add_argument('--orthodb-annotation-file', metavar='file', 
        help='File in uniprot_accession\\t format', type=argparse.FileType('r') )
    parser.add_argument('--orthodb-annotation-date', help='Date for the orthodb annotations')
    parser.add_argument('--annotation-output-file', metavar='file', 
        help='File that can be  copied into the db', type=argparse.FileType('w') )
    parser.add_argument('--annotation-source', help='Strongly recommended designation of the source for annotations', default='unknown-source')
    args = parser.parse_args()

    if args.obo_xml_file:
        data = parse_obo_xml(args.obo_xml_file)
        load_ontology(name=data['name'], terms=data['terms'], 
        subset_membership=data['subset_membership'], xrefs=data['xrefs'], 
        subsets=data['subsets'])

    if args.go_annotation_file:
        data = parse_go_annotation_file(args.go_annotation_file)
        if args.annotation_output_file:
            create_annotations(data, args.annotation_source, args.annotation_output_file)
        else:
            load_annotations(data, args.annotation_source)

    if args.orthodb_annotation_file:
        data = parse_orthodb_annotation_file(args.orthodb_annotation_file, args.orthodb_annotation_date)
        if args.annotation_output_file:
            create_annotations(data, args.annotation_source, args.annotation_output_file)
