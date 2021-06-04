# holds annotation objects
# e.g. TF family, expression, GO data


class Annotation():
    """Represents a flexible class that holds annotation info."""
    pass


# TODO: test this class
class Expression():
    def __init__(self, expr_dict):
        self.expr_dict = expr_dict # tissue -> rpkm
        self.avg_expr = self.compute_avg_expr()

    def __getitem__(self, tiss):
        # when expr_obj fetch by key (tissue), return value
        return self.expr_dict[tiss]


    def compute_avg_expr(self):
        tot = 0
        for k, v in self.expr_dict.items():
            v = float(v)
            tot += v
        avg_expr = tot/len(self.expr_dict)
        return avg_expr
