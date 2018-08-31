from circuit_block import Blocks

ps, pm, pg, eta, a0, a1, theta = \
(0.003, 0.003, 0.003, 0.01, 1.5, 0.0125, 0.24)

test_block = Blocks(ps, pm, pg, eta, a0, a1, theta)

p_succ = 0.1

n_attempts_lst = [test_block._success_number_of_attempts(p_succ) 
                    for _ in range(10**4)]
