import modules.vvMixinCore as vvMixinCore

class Validation():
    '''
    #Validation objects contain a number of variant interpretations
    '''
    pass

class ValOutput():
    '''
    #This object contains a single possible interpretation of a variant
    '''
    pass

class Validator(vvMixinCore.Mixin):
    '''
    #Mixins are used to split this very large, complex object over multiple files.
    #There is a logical chain to it, though:
    # vvMixinInit
    #     v
    # vvMixinConverters
    #     v
    # vvMixinCore
    #     v
    # Validator    <- this object.
    '''
    pass

