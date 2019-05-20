from VariantValidator import update_vv_db
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--delete', '-d', action='store_true', help='Delete the contents of the current database '
                                                                    'before updating')

    args = parser.parse_args()
    if args.delete:
        print("Deleting current database contents")
        update_vv_db.delete()

    update_vv_db.update()
