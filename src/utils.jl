"""
    @argcheck cond text

Throw an `ArgumentError` if `cond` is `false`.

adapted from base Julia's @assert macro.
"""
macro argcheck(ex, msgs...)
    msg = isempty(msgs) ? ex : msgs[1]
    if isa(msg, AbstractString)
        msg = msg # pass-through
    elseif !isempty(msgs) && (isa(msg, Expr) || isa(msg, Symbol))
        # message is an expression needing evaluating
        msg = :(Main.Base.string($(esc(msg))))
    elseif isdefined(Main, :Base) && isdefined(Main.Base, :string) && applicable(Main.Base.string, msg)
        msg = Main.Base.string(msg)
    else
        # string() might not be defined during bootstrap
        msg = quote
            msg = $(Expr(:quote,msg))
            isdefined(Main, :Base) ? Main.Base.string(msg) :
                (Core.println(msg); "Error during bootstrap. See stdout.")
        end
    end
    return :($(esc(ex)) ? $(nothing) : throw(ArgumentError($msg)))
end
