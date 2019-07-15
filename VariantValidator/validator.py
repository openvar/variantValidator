from .modules import vvMixinCore as vvMixinCore


class Validator(vvMixinCore.Mixin):
    """
    #Mixins are used to split this very large, complex object over multiple files.
    #There is a logical chain to it, though:
    # vvMixinInit
    #     v
    # vvMixinConverters
    #     v
    # vvMixinCore
    #     v
    # Validator    <- this object.
    """
    pass



