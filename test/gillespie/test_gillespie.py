from app.gillespie import gillespie
from unittest import TestCase

class TestGillespie(TestCase):
    def test_get_next_time_reaction(self):
        v = [1, 1]
        gillespie_res = gillespie.get_next_time_reaction(v)
        self.assertIsNotNone(gillespie_res.next_time)