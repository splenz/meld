import unittest
import mock


from meld import system
from meld.system import restraints


class TestAlwaysActiveCollection(unittest.TestCase):
    def setUp(self):
        self.coll = restraints.AlwaysActiveCollection()

    def test_adding_non_restraint_raises(self):
        with self.assertRaises(RuntimeError):
            self.coll.add_restraint(object)

    def test_should_be_able_to_add_restraint(self):
        rest = restraints.Restraint()
        self.coll.add_restraint(rest)
        self.assertIn(rest, self.coll.restraints)


class TestSelectivelyActiveCollection(unittest.TestCase):
    def test_restraint_should_be_present_after_adding(self):
        rest = [restraints.SelectableRestraint()]
        coll = restraints.SelectivelyActiveCollection(rest, 1)
        self.assertEqual(len(coll.restraints), 1)

    def test_can_add_two_restraints(self):
        rest = [restraints.SelectableRestraint(), restraints.SelectableRestraint()]
        coll = restraints.SelectivelyActiveCollection(rest, 1)
        self.assertEqual(len(coll.restraints), 2)

    def test_adding_non_selectable_restraint_should_raise(self):
        rest = [restraints.NonSelectableRestraint()]
        with self.assertRaises(RuntimeError):
            restraints.SelectivelyActiveCollection(rest, 1)

    def test_empty_restraint_list_should_raise(self):
        with self.assertRaises(RuntimeError):
            restraints.SelectivelyActiveCollection([], 0)

    def test_negative_num_active_should_raise(self):
        rest = [restraints.SelectableRestraint()]
        with self.assertRaises(RuntimeError):
            restraints.SelectivelyActiveCollection(rest, -1)

    def test_num_active_greater_than_num_restraints_should_raise(self):
        rest = [restraints.SelectableRestraint()]
        with self.assertRaises(RuntimeError):
            restraints.SelectivelyActiveCollection(rest, 2)

    def test_num_active_should_be_set(self):
        rest = [restraints.SelectableRestraint()]
        coll = restraints.SelectivelyActiveCollection(rest, 1)
        self.assertEqual(coll.num_active, 1)

    def test_should_wrap_bare_restraint_in_group(self):
        rest = [restraints.SelectableRestraint()]
        with mock.patch('meld.system.restraints.RestraintGroup.__init__', spec=True) as group_init:
            group_init.return_value = None
            restraints.SelectivelyActiveCollection(rest, 1)
            self.assertEqual(group_init.call_count, 1)

    def test_should_not_wrap_a_group_in_a_group(self):
        rest = [restraints.SelectableRestraint()]
        grps = [restraints.RestraintGroup(rest, 1)]
        with mock.patch('meld.system.restraints.RestraintGroup.__init__', spec=True) as group_init:
            restraints.SelectivelyActiveCollection(grps, 1)
            self.assertEqual(group_init.call_count, 0)


class TestRestraintGroup(unittest.TestCase):
    def test_should_accept_selectable_restraint(self):
        rest = [restraints.SelectableRestraint()]
        grp = restraints.RestraintGroup(rest, 1)
        self.assertEqual(len(grp.restraints), 1)

    def test_should_not_accept_non_selectable_restraint(self):
        rest = [restraints.NonSelectableRestraint()]
        with self.assertRaises(RuntimeError):
            restraints.RestraintGroup(rest, 1)

    def test_should_raise_on_empy_restraint_list(self):
        with self.assertRaises(RuntimeError):
            restraints.RestraintGroup([], 0)

    def test_num_active_below_zero_should_raise(self):
        rest = [restraints.SelectableRestraint()]
        with self.assertRaises(RuntimeError):
            restraints.RestraintGroup(rest, -1)

    def test_num_active_above_n_rest_should_raise(self):
        rest = [restraints.SelectableRestraint()]
        with self.assertRaises(RuntimeError):
            restraints.RestraintGroup(rest, 2)

    def test_num_active_should_be_set(self):
        rest = [restraints.SelectableRestraint()]
        grp = restraints.RestraintGroup(rest, 1)
        self.assertEqual(grp.num_active, 1)


class TestRestraintManager(unittest.TestCase):
    def setUp(self):
        self.mock_system = mock.Mock(spec=system.System)
        self.mock_system.index_of_atom.return_value = 0
        self.rest_manager = restraints.RestraintManager(self.mock_system)

    def test_can_add_as_always_active_non_selectable_restraint(self):
        rest = restraints.NonSelectableRestraint()
        self.rest_manager.add_as_always_active(rest)
        self.assertIn(rest, self.rest_manager.always_active)

    def test_can_add_as_always_active_selectable_restraint(self):
        rest = restraints.SelectableRestraint()
        self.rest_manager.add_as_always_active(rest)
        self.assertIn(rest, self.rest_manager.always_active)

    def test_can_add_list_of_always_active_restraints(self):
        rests = [restraints.SelectableRestraint(), restraints.NonSelectableRestraint()]
        self.rest_manager.add_as_always_active_list(rests)
        self.assertEqual(len(self.rest_manager.always_active), 2)

    def test_creating_bad_restraint_raises_error(self):
        with self.assertRaises(RuntimeError):
            self.rest_manager.create_restraint('blarg', x=42, y=99, z=-403)

    def test_can_create_distance_restraint(self):
        rest = self.rest_manager.create_restraint(
            'distance', atom_1_res_index=1, atom_1_name='CA',
            atom_2_res_index=2, atom_2_name='CA',
            r1=0, r2=0, r3=0.3, r4=999., k=2500)
        self.assertTrue(isinstance(rest, restraints.DistanceRestraint))

    def test_can_add_seletively_active_collection(self):
        rest_list = [restraints.SelectableRestraint(), restraints.SelectableRestraint()]
        self.rest_manager.add_selectively_active_collection(rest_list, 2)
        self.assertEqual(len(self.rest_manager.selectively_active_collections), 1)

    def test_can_create_restraint_group(self):
        rest_list = [restraints.SelectableRestraint(), restraints.SelectableRestraint()]
        grp = self.rest_manager.create_restraint_group(rest_list, 2)
        self.assertEqual(len(grp.restraints), 2)


class TestDistanceRestraint(unittest.TestCase):
    def setUp(self):
        self.mock_system = mock.Mock()

    def test_should_find_two_indices(self):
        restraints.DistanceRestraint(self.mock_system, 1, 'CA', 2, 'CA', 0, 0, 0.3, 999., 1.0)
        calls = [
            mock.call(1, 'CA'),
            mock.call(2, 'CA')]
        self.mock_system.index_of_atom.assert_has_calls(calls)

    def test_should_raise_on_bad_index(self):
        self.mock_system.index_of_atom.side_effect = KeyError()

        with self.assertRaises(KeyError):
            restraints.DistanceRestraint(self.mock_system, 1, 'BAD', 2, 'CA', 0, 0, 0.3, 999., 1.0)

    def test_should_raise_if_r2_less_than_r1(self):
        with self.assertRaises(RuntimeError):
            restraints.DistanceRestraint(self.mock_system, 1, 'CA', 2, 'CA', 10., 0., 10., 10., 1.0)

    def test_should_raise_if_r3_less_than_r2(self):
        with self.assertRaises(RuntimeError):
            restraints.DistanceRestraint(self.mock_system, 1, 'CA', 2, 'CA', 10., 10., 0., 10., 1.0)

    def test_should_raise_if_r4_less_than_r3(self):
        with self.assertRaises(RuntimeError):
            restraints.DistanceRestraint(self.mock_system, 1, 'CA', 2, 'CA', 10., 10., 10., 0., 1.0)

    def test_should_raise_with_negative_r(self):
        with self.assertRaises(RuntimeError):
            restraints.DistanceRestraint(self.mock_system, 1, 'CA', 2, 'CA', -1., 10., 10., 10., 1.0)

    def test_should_raise_with_negative_k(self):
        with self.assertRaises(RuntimeError):
            restraints.DistanceRestraint(self.mock_system, 1, 'CA', 2, 'CA', 10., 10., 10., 10., -1.0)


class TestTorsionRestraint(unittest.TestCase):
    def setUp(self):
        self.mock_system = mock.Mock()
        self.mock_system.index_of_atom.side_effect = [0, 1, 2, 3]

    def test_should_find_four_indices(self):
        restraints.TorsionRestraint(
            self.mock_system,
            1, 'CA',
            2, 'CA',
            3, 'CA',
            4, 'CA',
            180., 0., 1.0)

        calls = [
            mock.call(1, 'CA'),
            mock.call(2, 'CA'),
            mock.call(3, 'CA'),
            mock.call(4, 'CA')]
        self.mock_system.index_of_atom.assert_has_calls(calls)

    def test_should_raise_with_non_unique_indices(self):
        self.mock_system.index_of_atom.side_effect = [0, 0, 1, 2]  # repeated index
        with self.assertRaises(RuntimeError):
            restraints.TorsionRestraint(
                self.mock_system,
                1, 'CA',
                1, 'CA',
                3, 'CA',
                4, 'CA',
                180., 0., 1.0)

    def test_should_fail_with_phi_below_minus_180(self):
        with self.assertRaises(RuntimeError):
            restraints.TorsionRestraint(
                self.mock_system,
                1, 'CA',
                2, 'CA',
                3, 'CA',
                4, 'CA',
                -270., 0., 1.0)

    def test_should_fail_with_phi_above_180(self):
        with self.assertRaises(RuntimeError):
            restraints.TorsionRestraint(
                self.mock_system,
                1, 'CA',
                2, 'CA',
                3, 'CA',
                4, 'CA',
                270., 0., 1.0)

    def test_should_fail_with_delta_phi_above_180(self):
        with self.assertRaises(RuntimeError):
            restraints.TorsionRestraint(
                self.mock_system,
                1, 'CA',
                2, 'CA',
                3, 'CA',
                4, 'CA',
                0., 200., 1.0)

    def test_should_fail_with_delta_phi_below_0(self):
        with self.assertRaises(RuntimeError):
            restraints.TorsionRestraint(
                self.mock_system,
                1, 'CA',
                2, 'CA',
                3, 'CA',
                4, 'CA',
                0., -90., 1.0)

    def test_should_fail_with_negative_k(self):
        with self.assertRaises(RuntimeError):
            restraints.TorsionRestraint(
                self.mock_system,
                1, 'CA',
                2, 'CA',
                3, 'CA',
                4, 'CA',
                0., 90., -1.0)

