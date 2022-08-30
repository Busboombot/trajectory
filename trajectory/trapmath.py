X_LOWER_LIMIT = 10


def maxmin(l, v, h):
    return max(min(v, h), l)





def same_sign(a, b):
    return int(a) == 0 or int(b) == 0 or sign(a) == sign(b)


def triangular_area(v_0, v_1, v_c, a_max):
    """Area of a triangular profile, which accelerates to v_c and
    immediately decelerates"""
    return (2 * a_max * (v_0 + v_1) + (v_0 - v_c) ** 2 + (v_1 - v_c) ** 2) / (2 * a_max)


#####
# Group operations
#####




