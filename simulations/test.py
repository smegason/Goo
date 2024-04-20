def a(x=1):
    print(x)


def b(**kwargs):
    args = {"x": 1}
    args.update(kwargs)
    a(**args)


b()
