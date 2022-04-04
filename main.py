import pint

from parameters import ureg
import state
import integrated_model


def main(time: pint.Quantity, state: dict):
    model = integrated_model.MBRModel(state)
    t = 0 * ureg.s
    t_step = 60 * 15 * ureg.s
    while t < time:
        days = t.to('day').magnitude
        print(f'\r{days:.3f} days                     ', end='')
        state = model.step_model(t_step)
        t += t_step


if __name__ == "__main__":
    main(7 * ureg.day, state.starting_state)
