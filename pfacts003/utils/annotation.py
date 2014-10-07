'''
    Some weird thing for trying to deal with consensus annotations
    
    This is probably horribly horribly wrong.  Please don't track me down and
    make me defend it.  I will probably try to feed you a cookie in order to 
    distract you,
'''

# Horrible fix for 2.4
# Yes, this is bad.  make wsgi work with 2.7
try:
    any
except NameError:
    def any(s):
        for v in s:
            if v:
                return True
        return False

class Annotation:
    def __init__(self, annotation):
        self.annotation = annotation

    @classmethod
    def consensus(cls, annotations):
        return reduce(cls.combine_two_annotations, annotations)

    @classmethod
    def summarize(cls, annotations):
        return sorted(reduce(cls.add_annotation, annotations, {}).values())
        
    @classmethod
    def add_annotation(cls, set, annotation):
        if set.has_key(annotation.key()):
            pass
        else:
            set[annotation.key()] = annotation
        return set

class TreeAnnotation(Annotation):
    def key(self):
        return str(self.annotation)
    # annotation is an array of dicts
    # the array is ordered with the vaguest thing first and the
    # most refined last
    @classmethod
    def combine_two_annotations(cls, a1, a2):
        # first iteration
        if not a1:
            return a2

        output = []
        for a1_part, a2_part in zip(a1.annotation, a2.annotation):
            if a1_part == a2_part:
                output.append(a1_part)
            else:
                break
        
        return TreeAnnotation(output)

'''
    Yes, I know a set isn't ordered.

    The annotations are unique, but they are displayed in a specific order

    The concrete example is that we only want the best evidence code for a GO annotation
'''
class OrderedSetAnnotation(Annotation):
    def key(self):
        raise #needs to be subclassed

    @classmethod
    def add_annotation(cls, set, annotation):
        if set.has_key(annotation.key()):
            if annotation < set[annotation.key()]:
                set[annotation.key()] = annotation
        else:
            set[annotation.key()] = annotation
        return set


class GOAnnotation(OrderedSetAnnotation):
    def key(self):
        return self.annotation["accession"]
    def __lt__(self, other):
        return self.annotation["evidence_priority"] < other.annotation["evidence_priority"]
    def __gt__(self, other):
        return self.annotation["evidence_priority"] > other.annotation["evidence_priority"]



# Here are some functions that make those things useful
def run_class_method_on_homogenous_annotations(annotations, method_name):
    if not len(annotations):
        return annotations

    name = annotations[0].__class__.__name__

    if any(a.__class__.__name__ != name for a in annotations):
        raise Exception("List has conflicting annotation types")

    return getattr(annotations[0].__class__, method_name)(annotations)

def consensus(annotations):
    return run_class_method_on_homogenous_annotations(annotations, 'consensus')    

def summarize(annotations):
    return run_class_method_on_homogenous_annotations(annotations, 'summarize')    

